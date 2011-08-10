# Biolinux tests
#
#

def test_for_completed_install(ssh, box_name)
  # test for completion of the fabfile
  run "scp default:/var/log/biolinux.log biolinux.log"
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
