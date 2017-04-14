curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda_py3.sh
bash miniconda_py3.sh -b -p $HOME/miniconda3
echo "export PATH=\"$HOME/miniconda3/bin:\$PATH\"" >>$HOME/.bashrc
source $HOME/.bashrc		
conda update conda -y
conda install python=3.5 -y
conda config --add channels bioconda
conda config --add channels r
conda config --add channels ccbbucsd
conda install dual_crispr -y
conda install bioconductor-qvalue -y
set_up_dual_crispr