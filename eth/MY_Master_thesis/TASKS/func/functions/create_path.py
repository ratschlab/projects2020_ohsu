import os

def create_path(path,array):
    new_path=path
    for directory in array:
        new_path = os.path.join(new_path,directory)
        if not os.path.exists(new_path) and directory!=array[::-1][0]:
            os.mkdir(new_path)
    return new_path