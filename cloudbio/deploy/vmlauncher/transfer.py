import os
import gzip

from operator import itemgetter
from sys import exit
from threading import Thread
from threading import Condition
from Queue import Queue

from fabric.api import local, put, sudo, cd
from fabric.colors import red



class FileSplitter:
    """
    Works like the UNIX split command break up a file into parts like:
        filename_aaaaaaaaa
        filename_aaaaaaaab
        etc...
    """

    def __init__(self, chunk_size, destination_directory, callback):
        self.chunk_size = chunk_size * 1024 * 1024
        self.destination_directory = destination_directory
        self.chunk_callback = callback

    def split_file(self, path, compress, transfer_target):
        basename = os.path.basename(path)
        file_size = os.path.getsize(path)
        total_bytes = 0
        chunk_num = 0
        suffix = ''
        if compress:
            suffix = '.gz'

        input = open(path, 'rb')
        while True:
            chunk_name = "%s_part%08d%s" % (basename, chunk_num, suffix)
            chunk_path = os.path.join(self.destination_directory, chunk_name)
            this_chunk_size = min(self.chunk_size, file_size - total_bytes)
            if this_chunk_size <= 0:
                break

            chunk = input.read(this_chunk_size)
            total_bytes += len(chunk)
            if compress:
                chunk_output = gzip.open(chunk_path, 'wb')
            else:
                chunk_output = file(chunk_path, 'wb')
            chunk_output.write(chunk)
            chunk_output.close()

            self.chunk_callback.handle_chunk(chunk_path, transfer_target)
            chunk_num += 1


class TransferTarget:

    def __init__(self, file, precompressed, transfer_manager):
        self.file = file
        self.precompressed = precompressed
        self.do_compress = transfer_manager.compress
        self.do_split = transfer_manager.chunk_size > 0
        self.local_temp = transfer_manager.local_temp
        basename = os.path.basename(file)
        if len(basename) < 1:
            print red(Exception("Invalid file specified - %s" % file))
            exit(-1)
        self.basename = basename

    def should_compress(self):
        return not self.precompressed and self.do_compress

    def split_up(self):
        return self.do_split

    def clean(self):
        if self.should_compress():
            local("rm -rf '%s'" % self.compressed_file())

    def compressed_basename(self):
        if not self.precompressed:
            compressed_basename = "%s.gz" % self.basename
        else:
            compressed_basename = self.basename
        return compressed_basename

    def decompressed_basename(self):
        basename = self.basename
        if basename.endswith(".gz"):
            decompressed_basename = basename[:-len(".gz")]
        else:
            decompressed_basename = basename
        return decompressed_basename

    def compressed_file(self):
        compressed_file = "%s/%s.gz" % (self.local_temp, self.basename)
        return compressed_file

    def build_simple_chunk(self):
        if self.should_compress():
            compressed_file = self.compressed_file()
            local("gzip -f -9 '%s' -c > '%s'" % (self.file, compressed_file))
            return TransferChunk(compressed_file, self)
        else:
            return TransferChunk(self.file, self)


class TransferChunk:

    def __init__(self, chunk_path, transfer_target):
        self.chunk_path = chunk_path
        self.transfer_target = transfer_target

    def clean_up(self):
        was_split = self.transfer_target.split_up()
        was_compressed = self.transfer_target.should_compress()
        if was_split or was_compressed:
            local("rm '%s'" % self.chunk_path)


class FileTransferManager:

    def __init__(self,
                 compress=True,
                 num_compress_threads=1,
                 num_transfer_threads=1,
                 num_decompress_threads=1,
                 chunk_size=0,
                 transfer_retries=3,
                 destination="/tmp",
                 transfer_as="root",
                 local_temp=None):
        self.compress = compress
        self.num_compress_threads = num_compress_threads
        self.num_transfer_threads = num_transfer_threads
        self.num_decompress_threads = num_decompress_threads
        self.chunk_size = chunk_size
        self.transfer_retries = transfer_retries
        self.destination = destination
        self.transfer_as = transfer_as
        self.local_temp = local_temp

        if not self.local_temp:
            self.local_temp = "/tmp"

        local("mkdir -p '%s'" % self.local_temp)
        self.file_splitter = FileSplitter(self.chunk_size, self.local_temp, self)

    def handle_chunk(self, chunk, transfer_target):
        self._enqueue_chunk(TransferChunk(chunk, transfer_target))

    def transfer_files(self, files=[], compressed_files=[]):
        self.transfer_complete = False
        self.transfer_complete_condition = Condition()

        self._setup_destination_directory()

        self._setup_workers()

        self._enqueue_files(files, compressed_files)

        self._wait_for_completion()

    def _setup_workers(self):
        self._setup_compress_threads()
        self._setup_transfer_threads()
        self._setup_decompress_threads()

    def _setup_destination_directory(self):
        sudo("mkdir -p %s" % self.destination)
        self._chown(self.destination)

    def _setup_compress_threads(self):
        self.compress_queue = Queue()
        self._launch_threads(self.num_compress_threads, self._compress_files)

    def _setup_decompress_threads(self):
        self.decompress_queue = Queue()
        self._launch_threads(self.num_decompress_threads, self._decompress_files)

    def _setup_transfer_threads(self):
        self.transfer_queue = Queue()  # For now just transfer one file at a time
        self._launch_threads(self.num_transfer_threads, self._put_files)

    def _launch_threads(self, num_threads, func):
        for thread_index in range(num_threads):
            t = Thread(target=func)
            t.daemon = True
            t.start()

    def _enqueue_files(self, files, compressed_files):
        transfer_targets = []

        for file in files:
            transfer_target = TransferTarget(file, False, self)
            transfer_targets.append(transfer_target)

        for compressed_file in compressed_files:
            transfer_target = TransferTarget(compressed_file, True, self)
            transfer_targets.append(transfer_target)

        transfer_targets = self._sort_transfer_targets(transfer_targets)
        for transfer_target in transfer_targets:
            self.compress_queue.put(transfer_target)

    def _sort_transfer_targets(self, transfer_targets):
        for i in range(len(transfer_targets)):
            transfer_target = transfer_targets[i]
            transfer_targets[i] = transfer_target, os.stat(transfer_target.file).st_size
        transfer_targets.sort(key=itemgetter(1), reverse=True)
        return  [transfer_target[0] for transfer_target in transfer_targets]

    def _wait_for_completion(self):
        self.compress_queue.join()
        self.transfer_queue.join()
        self.transfer_complete_condition.acquire()
        self.transfer_complete = True
        self.transfer_complete_condition.notifyAll()
        self.transfer_complete_condition.release()
        self.decompress_queue.join()

    def _compress_files(self):
        while True:
            try:
                transfer_target = self.compress_queue.get()
                file = transfer_target.file
                if self.chunk_size > 0:
                    should_compress = transfer_target.should_compress()
                    self.file_splitter.split_file(file, should_compress, transfer_target)
                    self.decompress_queue.put(transfer_target)
                else:
                    simple_chunk = transfer_target.build_simple_chunk()
                    self._enqueue_chunk(simple_chunk)
            except Exception as e:
                print red("Failed to compress a file to transfer")
                print red(e)
            finally:
                self.compress_queue.task_done()

    def _decompress_files(self):
        if self.chunk_size > 0:
            self.transfer_complete_condition.acquire()
            while not self.transfer_complete:
                self.transfer_complete_condition.wait()
            self.transfer_complete_condition.release()
        while True:
            try:
                transfer_target = self.decompress_queue.get()
                basename = transfer_target.basename
                chunked = transfer_target.split_up()
                compressed = transfer_target.do_compress or transfer_target.precompressed
                with cd(self.destination):
                    if compressed and chunked:
                        destination = transfer_target.decompressed_basename()
                        if transfer_target.precompressed:
                            sudo("cat '%s_part'* | gunzip -c > %s" % (basename, destination), user=self.transfer_as)
                        else:
                            sudo("zcat '%s_part'* > %s" % (basename, destination), user=self.transfer_as)
                        sudo("rm '%s_part'*" % (basename), user=self.transfer_as)
                    elif compressed:
                        sudo("gunzip -f '%s'" % transfer_target.compressed_basename(), user=self.transfer_as)
                    elif chunked:
                        sudo("cat '%s'_part* > '%s'" % (basename, basename), user=self.transfer_as)
                        sudo("rm '%s_part'*" % (basename), user=self.transfer_as)
            except Exception as e:
                print red("Failed to decompress or unsplit a transfered file.")
                print red(e)
            finally:
                self.decompress_queue.task_done()

    def _put_files(self):
        while True:
            try:
                transfer_chunk = self.transfer_queue.get()
                transfer_target = transfer_chunk.transfer_target
                compressed_file = transfer_chunk.chunk_path
                basename = os.path.basename(compressed_file)
                self._put_as_user(compressed_file, "%s/%s" % (self.destination, basename))
                if not transfer_target.split_up():
                    self.decompress_queue.put(transfer_target)
            except Exception as e:
                print red("Failed to upload a file.")
                print red(e)
            finally:
                transfer_chunk.clean_up()
                self.transfer_queue.task_done()

    def _chown(self, destination):
        sudo("chown %s:%s '%s'" % (self.transfer_as, self.transfer_as, destination))

    def _put_as_user(self, source, destination):
        for attempt in range(self.transfer_retries):
            retry = False
            try:
                put(source, destination, use_sudo=True)
                self._chown(destination)
            except BaseException as e:
                retry = True
                print red(e)
                print red("Failed to upload %s on attempt %d" % (source, attempt + 1))
            except:
                # Should never get here, delete this block when more confident
                retry = True
                print red("Failed to upload %s on attempt %d" % (source, attempt + 1))
            finally:
                if not retry:
                    return
        print red("Failed to transfer file %s, exiting..." % source)
        exit(-1)

    def _enqueue_chunk(self, transfer_chunk):
        self.transfer_queue.put(transfer_chunk)
