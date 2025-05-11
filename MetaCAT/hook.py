import os
from stat import S_IRGRP, S_IROTH, S_IRUSR, S_IWUSR, S_IXGRP, S_IXOTH, S_IXUSR

from setuptools.command.build_py import build_py

flag = S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH


class PermissionHook(build_py):
    def run(self):
        super().run()
        for i in (
            'FragGeneScan-Darwin-arm64', 'FragGeneScan-Linux-x86_64',
            'hmmpress-Darwin-arm64', 'hmmpress-Linux-x86_64',
            'hmmsearch-Darwin-arm64', 'hmmsearch-Linux-x86_64',
            'mash-Darwin-arm64', 'mash-Linux-x86_64'
        ):
            path = os.path.join(self.build_lib, 'MetaCAT', i)
            os.chmod(path, flag)
        return None
