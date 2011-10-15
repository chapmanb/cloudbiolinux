# Biolinux shared tests

def test_for_completed_install(ssh, box_name)
  # test for completion of the fabfile
  run "scp vagrant:/var/log/biolinux.log biolinux.log"
  lastline = ssh.exec!("tail -1 /var/log/biolinux.log")
  print lastline
  error "Installation looks incomplete "+lastline if lastline !~ /#{box_name}/
  true
end

# Simple comparison using a regex
def test_for_match(ssh, cmd, regex)
  res = ssh.exec!(cmd)
  if res !~ /#{regex}/
    error(regex+" ("+cmd+") does not match "+res)
  end
end

def test_in_path(ssh, program)
  test_for_match(ssh,"which "+program, program)
end

def test_phylogeny_flavor(vagrant)
  vagrant.primary_vm.ssh.execute do |ssh|
    # test for command line tools
    test_for_match(ssh,"beast -help", "BEAST")
    test_in_path(ssh,"mb")
    test_in_path(ssh,"mpirun.openmpi")
    test_in_path(ssh,"mrbayes-multi")
    # test_for_match(ssh,"mb --help", "MrBayes")
    # test_for_match(ssh,"mpirun.openmpi -np 4 mrbayes-multi -h", "MrBayes")
    test_for_match(ssh,"raxmlHPC -v", "RAxML")
    test_in_path(ssh,"codeml")
  end 
end
