# standard libraries
import argparse
import os
import shutil

from dual_crispr import dual_crispr_pipeliner as ns_dcpipe


def _parse_cmd_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="path to config file; if not specified, config file in default location will "
                                         "be used")
    args = parser.parse_args()
    return args.config


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
    config_fp = _parse_cmd_line_args()
    params = ns_dcpipe.set_up_data_subdirs_and_get_machine_configs(config_fp)
    source_dir = os.path.join(params["main_dir"], "test_data")
    dest_dir = params[ns_dcpipe.DirectoryKeys.RAW_DATA.value]
    _copy_test_files(source_dir, dest_dir)


if __name__ == '__main__':
    main()