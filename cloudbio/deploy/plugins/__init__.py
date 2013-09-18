"""
Deploy plugins.

local_actions:
  These actions occur before a VM has been created, if only local actions are specified, no VM is launched.

configure_actions:
  Take in only options and use fabric to configure instance in some way.

ready_actions:
  Once VM is configured, these actions are executed.

compound_actions:
  Dictionary of list values. Key represents an short-cut action that is expanded to each action specified
  in corresponding value (a list of simple actions - standard or defined in the plugins).
"""
