import argparse
import os

from epidemiology.epidata import defaultDict as dd

def clean_data(all, rki, john_hopkins, spain, population, hdf5, out_path):

    ending = "json"
    if(hdf5):
        ending = "h5"

    if all:

        # TODO: make general dictionary with all countries used
        directories = ['Germany/', 'Spain/', 'France/', 'Italy/', 'US/', 'SouthKorea/', 'China/']

        # delete files in directory
        for cdir in directories:
            directory = os.path.join(out_path, cdir)

            try:
                files = os.listdir(directory)
            except FileNotFoundError:
                continue

            for item in files:
                print("Deleting file ", os.path.join(directory, item))
                os.remove(os.path.join(directory, item))

            # delete directories
            print("Deleting directory ", directory)
            os.rmdir(directory)

        # delete further jh files
        files = []
        try:
            files = os.listdir(out_path)
        except FileNotFoundError:
            pass

        for item in files:
            if item.endswith(".json") or item.endswith(".h5"):
                print("Deleting file ", os.path.join(out_path, item))
                os.remove(os.path.join(out_path, item))

    elif(rki):
        directory = os.path.join(out_path, 'Germany/')
        files = []
        try:
            files = os.listdir(directory)
        except FileNotFoundError:
            pass

        for item in files:
            if item.endswith(ending):
                if "_rki" in item or "RKI" in item:
                    print("Deleting file ", os.path.join(directory, item))
                    os.remove(os.path.join(directory, item))

        # delete directory if empty
        try:
            os.rmdir(directory)
        except OSError:
            return

        print("Deleting directory ", directory)

    elif (population):
        directory = os.path.join(out_path, 'Germany/')
        files = []
        try:
            files = os.listdir(directory)
        except FileNotFoundError:
            pass

        for item in files:
            if item.endswith(ending):
                if "Popul" in item or "FullDataB" in item or "FullDataL" in item:
                    print("Deleting file ", os.path.join(directory, item))
                    os.remove(os.path.join(directory, item))

        # delete directory if empty
        try:
            os.rmdir(directory)
        except OSError:
            return

        print("Deleting directory ", directory)


    elif (spain):
        directory = os.path.join(out_path, 'Spain/')
        files = []
        try:
            files = os.listdir(directory)
        except FileNotFoundError:
            pass

        for item in files:
            if item.endswith(ending) and "spain" in item:
                print("Deleting file ", os.path.join(directory, item))
                os.remove(os.path.join(directory, item))

        # delete directory if empty
        try:
            os.rmdir(directory)
        except OSError:
            return;

        print("Deleting directory ", directory)

    elif(john_hopkins):
        # TODO: make general dictionary with all countries used
        directories = ['Germany/', 'Spain/', 'France/', 'Italy/', 'US/', 'SouthKorea/', 'China/']

        # delete files in directory
        for cdir in directories:
            directory = os.path.join(out_path, cdir)

            try:
                files = os.listdir(directory)
            except FileNotFoundError:
                continue

            for item in files:
                if item.endswith(ending) and "_jh" in item:
                    print("Deleting file ", os.path.join(directory, item))
                    os.remove(os.path.join(directory, item))

            # delete directories
            try:
               os.rmdir(directory)
            except OSError:
                continue;

            print("Deleting directory ", directory)

        # delete further jh files
        files = []
        try:
            files = os.listdir(out_path)
        except FileNotFoundError:
            pass

        for item in files:
            if item.endswith(ending):
                if "_jh" in item or "JohnHopkins" in item:
                    print("Deleting file ", os.path.join(out_path, item))
                    os.remove(os.path.join(out_path, item))

    else:
        print("Please specify what should be deleted. See --help for details.")


def cli():

    out_path_default = dd.defaultDict['out_folder']

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--all',
                        help='Deletes all data and folders which could be possibly written by the epidata package.',
                        action='store_true')

    parser.add_argument('-r', '--rki', help='Deletes just rki data.',
                        action='store_true')
    parser.add_argument('-j', '--john-hopkins', help='Deletes just data from John Hopkins university.',
                        action='store_true')
    parser.add_argument('-s', '--spain', help='Deletes just spain data.',
                        action='store_true')
    parser.add_argument('-p', '--population', help='Deletes just population data.',
                        action='store_true')

    parser.add_argument('-h5', '--hdf5', help='Deletes just hdf5 files.',
                        action='store_true')
    parser.add_argument('-o', '--out_path', type=str, default=out_path_default, help='Defines folder for output.')

    args = parser.parse_args()

    return [args.all, args.rki, args.john_hopkins, args.spain, args.population, args.hdf5, args.out_path]

def main():

   [all, rki, john_hopkins, spain, population, hdf5, out_path] = cli()

   clean_data(all, rki, john_hopkins, spain, population, hdf5, out_path)

if __name__ == "__main__":

   main()