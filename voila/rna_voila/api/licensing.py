from cryptography.fernet import Fernet
from datetime import datetime
import sys, os

LICENSE_CHECKING_ENABLED=False
key = b'7ulXO1YfJWoPscAsklVwBtaHJKzfGQBlPS7Rj4DCoD8='

def _check_license(filepath):
    with open(filepath, 'rb') as f:
        data = f.read()

    cipher_suite = Fernet(key)

    try:
        decoded_text = cipher_suite.decrypt(data).decode()
    except:
        print(f"There was an error reading the license file {filepath}")
        sys.exit(1)

    label, datestr = decoded_text.split('_')
    if datestr == "None":
        return True
    try:
        exp_date = datetime.strptime(datestr, r'%d-%m-%Y')
    except:
        print(f"There was an error reading the license file {filepath}")
        sys.exit(1)
    cur_date = datetime.now()
    if cur_date > exp_date:
        print(f"The License {label} ({filepath}) (expires {exp_date.strftime(r'%d-%m-%Y')}) has expired, exiting")
        sys.exit(1)

    print(f"Using License {label} ({filepath}) (expires {exp_date.strftime(r'%d-%m-%Y')})")

    return True
def check_license(license_file_path):

    if not LICENSE_CHECKING_ENABLED:
        return True

    if license_file_path:
        _check_license(license_file_path)
    else:
        # check envvar
        license_file_envvar = os.getenv('MAJIQ_LICENSE_FILE')
        if license_file_envvar and os.path.exists(license_file_envvar):
            _check_license(license_file_envvar)
        else:
            # check local directory
            files = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('majiq_license')]
            if files:
                _check_license(files[0])
            else:
                # check home directory
                basepath = os.path.expanduser('~')
                paths = [os.path.join(basepath, f) for f in os.listdir(basepath)]
                files = [f for f in paths if os.path.isfile(f) and os.path.basename(f).startswith('majiq_license')]

                if files:
                    _check_license(files[0])
                else:
                    print(f'No Majiq License files were found, exiting')
                    sys.exit(1)



