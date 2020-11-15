
#!/usr/bin/env bash

USERNAME=$1
IP_ADDRESS=$2

if [[ -z "$USERNAME" ]]; then 
  echo "✋ missing username"
  echo "Usage: provision.sh <user name> <ip address>"
  exit
fi

if [[ -z "$IP_ADDRESS" ]]; then 
  echo "✋ missing IP"
  echo "Usage: provision.sh <user name> <ip address>"
  exit
fi

echo "🔨 Adding 👤 $USERNAME 💻 $IP_ADDRESS"
echo "sudo /home/ubuntu/projman/add_user.sh $USERNAME $IP_ADDRESS" | ssh  -i keys/projman_key ubuntu@$(terraform output | awk -F' = ' '{print $2}')
