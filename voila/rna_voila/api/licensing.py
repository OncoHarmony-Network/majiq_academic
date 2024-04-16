from cryptography.fernet import Fernet
from datetime import datetime
import sys, os

LICENSE_CHECKING_TYPE = "academic"
key = b'7ulXO1YfJWoPscAsklVwBtaHJKzfGQBlPS7Rj4DCoD8='


def td_format(td_object):
    seconds = int(td_object.total_seconds())
    periods = [
        ('year',  60*60*24*365),
        ('month', 60*60*24*30),
        ('day',   60*60*24),
    ]
    strings = []
    for period_name, period_seconds in periods:
        if seconds > period_seconds:
            period_value , seconds = divmod(seconds, period_seconds)
            has_s = 's' if period_value > 1 else ''
            strings.append("%s %s%s" % (period_value, period_name, has_s))

    return ", ".join(strings)

def _check_license(filepath, log_obj):
    with open(filepath, 'rb') as f:
        data = f.read()

    cipher_suite = Fernet(key)

    try:
        decoded_text = cipher_suite.decrypt(data).decode()
    except:
        print(f"There was an error reading the license file {filepath}")
        sys.exit(1)

    usage, label, datestr = decoded_text.split('_')
    if usage != LICENSE_CHECKING_TYPE:
        log_obj.critical(f"License file type '{usage}' doesn't match Majiq distribution '{LICENSE_CHECKING_TYPE}', exiting")
        sys.exit(1)

    if datestr == "None":
        exp_str = 'Never'
    else:
        try:
            exp_date = datetime.strptime(datestr, r'%d-%m-%Y')
        except:
            log_obj.critical(f"There was an error reading the license file {filepath}")
            sys.exit(1)
        cur_date = datetime.now()
        if cur_date > exp_date:
            log_obj.critical(f"The License {label} ({filepath}) (expires {exp_date.strftime(r'%d-%m-%Y')}) has expired, exiting")
            sys.exit(1)
        exp_str = f"{exp_date.strftime(r'%d-%m-%Y')} ({td_format(exp_date - cur_date)})"

    jst = 64
    log_obj.info("╔═══════════════════════════════════════════════════════════════╗")
    log_obj.info(f"╠╡ {usage.upper()} License applied".ljust(jst)[:jst] + "║")
    log_obj.info(f"║  Name: {label}".ljust(jst)[:jst] + "║")
    log_obj.info(f"║  File: {os.path.basename(filepath)}".ljust(jst)[:jst] + "║")
    log_obj.info(f"║  Expiration Date: {exp_str}".ljust(jst)[:jst] + "║")
    if usage == 'academic':
        log_obj.info("║  ".ljust(jst)[:jst] + "║")
        log_obj.info("╠╡ The academic license is for non-commercial purposes by  ".ljust(jst)[:jst] + "║")
        log_obj.info("╠╡ individuals at an academic or not for profit institution. ".ljust(jst)[:jst] + "║")
    log_obj.info("╚═══════════════════════════════════════════════════════════════╝")


    return True
def check_license(license_file_path, log_obj):

    if LICENSE_CHECKING_TYPE == "None":
        return True

    if license_file_path:
        _check_license(license_file_path, log_obj)
    else:
        # check envvar
        license_file_envvar = os.getenv('MAJIQ_LICENSE_FILE')
        if license_file_envvar and os.path.exists(license_file_envvar):
            _check_license(license_file_envvar, log_obj)
        else:
            # check local directory
            files = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('majiq_license')]
            if files:
                _check_license(files[0], log_obj)
            else:
                # check home directory
                basepath = os.path.expanduser('~')
                paths = [os.path.join(basepath, f) for f in os.listdir(basepath)]
                files = [f for f in paths if os.path.isfile(f) and os.path.basename(f).startswith('majiq_license')]

                if files:
                    _check_license(files[0], log_obj)
                else:
                    log_obj.critical(f'No Majiq License files were found, exiting')
                    sys.exit(1)



