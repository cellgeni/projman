
#!/usr/bin/env bash

USERNAME=$1
IP_ADDRESS=$2

if [[ -z "$USERNAME" ]] || [[ -z "$IP_ADDRESS" ]]; then 
  echo "✋ missing username or IP"
  echo "Usage: provision.sh <user name> <ip address>"
  exit
fi


echo "🔨 Adding 👤 $USERNAME 💻 $IP_ADDRESS"
echo "sudo /home/ubuntu/projman/add_user.sh $USERNAME $IP_ADDRESS" | ssh -o StrictHostKeyChecking=no -i keys/projman_key ubuntu@$(terraform output | awk -F' = ' '{print $2}')
