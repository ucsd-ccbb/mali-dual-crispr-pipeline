# standard libraries
import os
import shutil
import sys

import dual_crispr  # do not remove this "unused" import--it IS used in call to sys.modules["dual_crispr"]
import ccbb_pyutils.files_and_paths as ns_file


# code from https://stackoverflow.com/questions/1994488/copy-file-or-directories-recursively-in-python#1994840
def _copy_test_files(root_src_dir, root_dst_dir):
    for src_dir, dirs, files in os.walk(root_src_dir):
        dst_dir = src_dir.replace(root_src_dir, root_dst_dir, 1)

        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)

        for file_ in files:
            src_file = os.path.join(src_dir, file_)
            dst_file = os.path.join(dst_dir, file_)

            if os.path.exists(dst_file):
                os.remove(dst_file)
            shutil.copy(src_file, dst_dir)


def main():
    home_dir_path = os.environ['HOME']

    # get the location of the installed files
    pkg_dir_path = os.path.dirname(sys.modules["dual_crispr"].__file__)
    files_dir_path = os.path.join(pkg_dir_path, "distributed_files")
    home_dual_crispr_path = os.path.join(home_dir_path, "dual_crispr")
    ns_file.verify_or_make_dir(home_dual_crispr_path)

    # move them to the user's home directory
    _copy_test_files(files_dir_path, home_dual_crispr_path)


if __name__ == '__main__':
    main()