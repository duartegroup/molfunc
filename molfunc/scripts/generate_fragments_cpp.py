"""
Generate the .cpp with all the fragments from the library (src/data/)
"""
import os
import sys

root_dir = sys.argv[1]


def print_fragment_lib():
    """Print the fragment library"""

    with open(f'{root_dir}/src/species/fragments.cpp', 'a') as cpp_file:

        print('\n\nnamespace molfunc{',
              '     FragmentLib::FragmentLib(){',
              '         fragments = {',
              '',
              sep='\n', file=cpp_file)

        data_dir = os.path.join(root_dir, 'src', 'species', 'data')

        for filename in os.listdir(data_dir):
            if not filename.endswith('.xyz'):
                continue

            aliases = []
            xyzs = []

            for i, line in enumerate(open(os.path.join(data_dir, filename), 'r')):
                if i == 0:
                    continue

                if i == 1:
                    aliases = [item.strip()
                               for item in line.split(' ')[-1].split(',')]
                    continue

                if len(line.split()) == 0:
                    break

                try:
                    symbol, x, y, z = line.split()
                    xyzs.append((symbol, x, y, z))
                except ValueError:
                    raise ValueError('Mallformatted xyz line:'
                                     f'{line}')

            print('         Fragment({', file=cpp_file)
            for (symbol, x, y, z) in xyzs:
                print(f'            Atom3D("{symbol}", {x}, {y}, {z}),', file=cpp_file)
            print('         },\n'
                  '         {'+''.join([f'"{item}",' for item in aliases])+'})\n'
                  '         ,',
                  file=cpp_file)

        print('     };\n'
              '  }\n'
              '}',
              file=cpp_file)

    return None


def fragment_lib_present():
    """Is the fragment library already present in the source file?"""

    for line in open(f'{root_dir}/src/species/fragments.cpp', 'r'):
        if 'FragmentLib::FragmentLib()' in line:
            return True

    return False


if __name__ == '__main__':

    if not fragment_lib_present():
        print_fragment_lib()
