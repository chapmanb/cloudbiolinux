# Run a command line command, without intercepting stdout and stderr
def run cmd
  print "--> ",cmd
  res = system(cmd)  # don't use backticks for fab!
  if res == nil
    error "Command failed!"
  end
  res
end

# Error handler
def error msg
  $stderr.print "\nFATAL ERROR: ",msg,"\n"
  exit 1
end


