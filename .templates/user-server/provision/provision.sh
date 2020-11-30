#!/usr/bin/env bash

if [[ "$EUID" -ne 0 ]]; then
  echo "✋ run this script as root"
  exit
fi

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

function wait_for {
  i=0
  tput sc
  while fuser $WAITING_FOR >/dev/null 2>&1 ; do
      case $(($i % 8)) in
          0 ) j="🌑" ;;
          1 ) j="🌒" ;;
          2 ) j="🌓" ;;
          3 ) j="🌔" ;;
          4 ) j="🌕" ;;
          5 ) j="🌖" ;;
          6 ) j="🌗" ;;
          7 ) j="🌘" ;;
      esac
      tput rc
      echo -en "\r[$j] waiting for lock on $WAITING_FOR"
      sleep 1
      ((i=i+1))
  done
  sleep 2
}

# update packages
apt-get update
WAITING_FOR=/var/lib/dpkg/lock
wait_for
WAITING_FOR=/var/lib/apt/lists/lock
wait_for
WAITING_FOR=/var/lib/dpkg/lock-frontend
wait_for
# install dependencies
apt-get install -y -qq nginx libcrack2 ipython sshfs

# fix nginx permissions
chown -R ubuntu:ubuntu /etc/nginx/
chown -R ubuntu:ubuntu /var/www/html/

# nameserver resolve
cat > /etc/netplan/60-resolve.yaml <<EOF
network:
    version: 2
    ethernets:
        ens3:
            dhcp4: true
            nameservers:
               search: [internal.sanger.ac.uk]
EOF
netplan generate
netplan apply


###############################################################################
# projman master init
#

serverURL="projman.cellgeni.sanger.ac.uk"
echo "$IP_ADDRESS" > /projman/hostIP
echo "$USERNAME" > /projman/userName

# copy stuff where it's supposed to be
chmod +x /projman/mount-farm && mv /projman/mount-farm /usr/local/bin/
chmod +x /projman/projMan && mv /projman/projMan /usr/local/bin/
sed "s@<!-- TARGET FOR IP-->@$IP_ADDRESS;@g;" /projman/index.html > /var/www/html/index.html
sed "s@server_name.*@server_name ${serverURL};@g; s@base_location.*@location /$USERNAME/ {@g;" /projman/server.conf > /etc/nginx/conf.d/server.conf

nginx -s reload

# create mount-farm endpoints
mkdir -p /nfs
mkdir -p /lustre
mkdir -p /warehouse
#
###############################################################################