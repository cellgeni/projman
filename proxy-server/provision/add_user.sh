#!/usr/bin/env bash

if [[ "$EUID" -ne 0 ]]; then
  echo "✋ run this script as root"
  exit
fi

USERNAME=$1
IP_ADDRESS=$2

if [[ -z "$USERNAME" ]]; then 
  echo "✋ missing username"
  echo "Usage: ./add_user.sh <user name> <ip address>"
  exit
fi

if [[ -z "$IP_ADDRESS" ]]; then 
  echo "✋ missing IP"
  echo "Usage: ./add_user.sh <user name> <ip address>"
  exit
fi

read -r -d '' TEMPLATE << EOM
    location /$USERNAME/ {
        proxy_pass                  http://$IP_ADDRESS/$USERNAME/;
        proxy_set_header            Host \$host;
        proxy_set_header            X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header            X-Real-IP \$remote_addr;
        proxy_cookie_domain         localhost projman.cellgeni.sanger.ac.uk;
        proxy_set_header            X-NginX-Proxy true;
        proxy_ssl_session_reuse     off;
        proxy_redirect              off;
        proxy_buffering             off;
        proxy_http_version          1.1;
        proxy_set_header            Upgrade \$http_upgrade;
        proxy_set_header            Connection "upgrade";
    }
EOM

TEMPLATE=`echo ${TEMPLATE} | tr '\n' "\\n"`

sed -i "s@#NEW LOCATION PLACEHOLDER#@${TEMPLATE}\\n\\n#NEW LOCATION PLACEHOLDER#@g" /etc/nginx/conf.d/server.conf

nginx -s reload