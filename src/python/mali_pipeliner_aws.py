# standard libraries
import os

# ccbb libraries
from dual_crispr_pipeliner import DirectoryKeys, parse_cmd_line_args

# project-specific libraries
from mali_pipeliner import run_mali_pipeline

__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def main():
    fastq_set_name, human_readable_name, library_name = parse_cmd_line_args()

    # *********************************************************
    # you'll probably decide on these values before beginning screen analysis, set them once, and then leave them alone
    # unless you change your project folder and/or hardware
    num_processors = 3
    main_dir = "/home/ec2-user/jupyter-genomics"
    data_dir = "/data"

    library_def_dir = os.path.join(main_dir, "library_definitions")

    # *********************************************************
    # you'll probably never change these, assuming you use the suggested file structure
    pipeline_name = "crispr"

    dirs_dict = {}
    dirs_dict[DirectoryKeys.CODE] = os.path.join(main_dir, "src", pipeline_name)
    dirs_dict[DirectoryKeys.NOTEBOOKS] = os.path.join(main_dir, DirectoryKeys.NOTEBOOKS.value, pipeline_name)
    dirs_dict[DirectoryKeys.LIBRARIES] = # os.path.join(data_dir, DirectoryKeys.RAW_DATA.value, "general")
    dirs_dict[DirectoryKeys.RAW_DATA] = os.path.join(data_dir, DirectoryKeys.RAW_DATA.value)
    dirs_dict[DirectoryKeys.INTERIM_DATA] = os.path.join(data_dir, DirectoryKeys.INTERIM_DATA.value, fastq_set_name)
    dirs_dict[DirectoryKeys.PROCESSED_DATA] = os.path.join(data_dir, DirectoryKeys.PROCESSED_DATA.value)

    # *********************************************************
    # DON'T change this unless you *really* know what you're doing :)
    run_mali_pipeline(fastq_set_name, human_readable_name, library_name, num_processors, dirs_dict)


if __name__ == '__main__':
    main()
