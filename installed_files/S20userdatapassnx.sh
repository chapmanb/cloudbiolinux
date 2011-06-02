#!/bin/bash

# Check for passwords sent in via the user-data box in the AWS
# console. This prevents ever needing to SSH into the instance
# to allow connections via FreeNX directly.

#Get password from user data set in the AWS console
#The user data must be a single line beginning 'password='. 
UD=`curl --connect-timeout 5 -s http://169.254.169.254/latest/user-data | head -n1 | grep -i '^password='`

if [[ -n "$UD" ]]; then
        echo "Setting ubuntu password from user-data"
        #update the password for the ubuntu user
        /usr/bin/autopasswd ubuntu ${UD#*=} >&/dev/null
        
        #force SSH to allow password logins
        sed -i 's/^PasswordAuthentication .*/PasswordAuthentication yes/' /etc/ssh/sshd_config
        /etc/init.d/ssh reload >&/dev/null
        
        #set-up FreeNX
        dpkg-reconfigure -pcritical freenx-server >&/dev/null
fi

