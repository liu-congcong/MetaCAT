import os
from stat import S_IRGRP, S_IROTH, S_IRUSR, S_IWUSR, S_IXGRP, S_IXOTH, S_IXUSR

from setuptools.command.build_py import build_py

flag = S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH


class PermissionHook(build_py):
    def run(self):
        super().run()
        path = os.path.join(self.build_lib, 'MetaCAT')
        dirEntry = os.scandir(path)
        for i in dirEntry:
            if i.is_file() and i.name.endswith(('darwin-arm64', 'linux-x86_64', 'windows-amd64', 'windows-x86_64')):
                os.chmod(i.path, flag)
        dirEntry.close()
        return None
