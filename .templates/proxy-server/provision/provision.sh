#!/usr/bin/env bash

if [[ "$EUID" -ne 0 ]]; then
  echo "✋ run this script as root"
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
apt-get install -y -qq nginx

# fix nginx permissions
chown -R ubuntu:ubuntu /etc/nginx/conf.d/

# move stuff
mv /home/ubuntu/projman/server.conf /etc/nginx/conf.d/server.conf

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

