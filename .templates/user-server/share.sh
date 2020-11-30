MACHINE_IP=$(terraform output | grep server_ip | awk -F' = ' '{print $2}')
USERNAME=$(terraform output | grep server_name | awk -F' = ' '{print $2}')

# copy keys and create readme
echo "🗃️ gathering files"
mkdir -p tmp
cp keys/user_key tmp/${USERNAME}.projman.key
cp keys/user_key.pub tmp/${USERNAME}.projman.key.pub
cat > tmp/README.txt <<EOF
# JupyterHub project launcher
To manage existing projects, or create new projects ssh into your machine at ${MACHINE_IP} and run "projMan".

# Connect using
ssh -i ${USERNAME}.projman.key ubuntu@${MACHINE_IP}
EOF

# tar files and cleanup
echo "📦 packing"
tar -cjf ${USERNAME}_projman.tar  -C tmp .
echo "🧹 cleanup"
rm -rf tmp

echo "📨 share ${USERNAME}_projman.tar !"