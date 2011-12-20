# Run a command line command, without intercepting stdout and stderr
def run cmd
  print "--> ",cmd
  system(cmd)  # don't use backticks for fab!
end

# Error handler
def error msg
  $stderr.print "FATAL ERROR: ",msg
  exit 1
end


