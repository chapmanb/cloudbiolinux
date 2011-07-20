# Biolinux tests
#
#

def test_for_completed_install(edition_name)
  # test for completion of the fabfile
  run "scp vagrant:/var/log/biolinux.log biolinux.log"
  lastline = `tail -1 biolinux.log`
  print lastline
  error "Installation looks incomplete "+lastline if lastline !~ /#{edition_name}/
  true
end

# Simple comparison using a regex
def test_for_match(ssh, cmd, regex)
  res = ssh.exec!(cmd)
  if res !~ /#{regex}/
    error(regex+" ("+cmd+") does not match "+res)
  end
end
