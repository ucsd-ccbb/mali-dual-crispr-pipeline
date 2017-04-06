# Jupyter set-up script based on original by Cathy Yeh, from
# http://efavdb.com/deep-learning-with-jupyter-on-aws/

CERTIFICATE_DIR="$HOME/certificate"
JUPYTER_CONFIG_DIR="$HOME/.jupyter"

if [ ! -d "$CERTIFICATE_DIR" ]; then
    mkdir $CERTIFICATE_DIR
    openssl req -x509 -nodes -days 365 -newkey rsa:1024 -keyout "$CERTIFICATE_DIR/mykey.key" -out "$CERTIFICATE_DIR/mycert.pem" -batch
    chown -R ec2-user $CERTIFICATE_DIR
fi

if [ ! -f "$JUPYTER_CONFIG_DIR/jupyter_notebook_config.py" ]; then
    mkdir $JUPYTER_CONFIG_DIR

    # append notebook server settings
    cat <<EOF >> "$JUPYTER_CONFIG_DIR/jupyter_notebook_config.py"
# Set options for certfile, ip, password, and toggle off browser auto-opening
c.NotebookApp.certfile = u'$CERTIFICATE_DIR/mycert.pem'
c.NotebookApp.keyfile = u'$CERTIFICATE_DIR/mykey.key'
# Set ip to '*' to bind on all interfaces (ips) for the public server
c.NotebookApp.ip = '*'
c.NotebookApp.password = u'sha1:a4d000f86fc5:b550b1686b71e50ac87d52a72500034c84f5e3f9'
c.NotebookApp.open_browser = False
# It is a good idea to set a known, fixed port for server access
c.NotebookApp.port = 8888
EOF
    chown -R ec2-user $JUPYTER_CONFIG_DIR
else
    echo "Error in dual_crispr run_set_up_pipeline: Cannot create jupyter config file; $JUPYTER_CONFIG_DIR/jupyter_notebook_config.py already exists."
fi