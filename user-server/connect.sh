ssh -i keys/user_key ubuntu@$(terraform output | grep server_ip | awk -F' = ' '{print $2}')
