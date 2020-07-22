
import os, tempfile, shutil, time
from subprocess import Popen, PIPE, STDOUT, check_output
from distutils.dir_util import copy_tree

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
import traceback

"""
!!!

To run this script::::

-you must have compiled and installed majiq such that the version you are testing is the one in the path that will
run if you run 'majiq'
-you must download 'chromedriver' for your OS/chrome version and place it in the same directory as this file and
it must be called 'chromedriver' (https://chromedriver.chromium.org/downloads)
-you must download the "workshop example data" (https://majiq.biociphers.org/download/workshop_example.zip) and put
the contents of it in the "data/workshop_data" folder adjacent to this file. 

!!!
"""

class TestWebDriver:

    def __init__(self):
        options = Options()
        options.headless = True

        self.driver = webdriver.Chrome(os.path.abspath(os.path.join(os.path.dirname(__file__), 'chromedriver')), chrome_options=options)

    def load_page(self, url):
        self.driver.get(url)
        while True:
            time.sleep(1)
            if self.page_has_loaded():
                break
        #time.sleep(13)

    def end(self):
        self.driver.close()

    def page_has_loaded(self):
        page_state = self.driver.execute_script('return document.readyState;')
        return page_state == 'complete'

    def get_num_records(self):
        records = int(self.driver.execute_script('return $results.page.info().recordsTotal;'))
        print(records)
        return records

    def get_performance_timing(self):
        return int(self.driver.execute_script('return window.performance.timing.domContentLoadedEventEnd - window.performance.timing.navigationStart;'))

    def get_all_records(self):
        return self.driver.execute_script('''
        await fetch("http://localhost:5006/index-table", {"credentials":"include","headers":{"accept":"application/json, text/javascript, */*; q=0.01","accept-language":"en-US,en;q=0.9","cache-control":"no-cache","content-type":"application/x-www-form-urlencoded; charset=UTF-8","pragma":"no-cache","x-requested-with":"XMLHttpRequest"},"referrer":"http://localhost:5006/","referrerPolicy":"no-referrer-when-downgrade","body":"draw=1&columns%5B0%5D%5Bdata%5D=0&columns%5B0%5D%5Bname%5D=&columns%5B0%5D%5Bsearchable%5D=true&columns%5B0%5D%5Borderable%5D=true&columns%5B0%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B0%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B1%5D%5Bdata%5D=1&columns%5B1%5D%5Bname%5D=&columns%5B1%5D%5Bsearchable%5D=true&columns%5B1%5D%5Borderable%5D=true&columns%5B1%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B1%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B2%5D%5Bdata%5D=2&columns%5B2%5D%5Bname%5D=&columns%5B2%5D%5Bsearchable%5D=true&columns%5B2%5D%5Borderable%5D=false&columns%5B2%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B2%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B3%5D%5Bdata%5D=3&columns%5B3%5D%5Bname%5D=&columns%5B3%5D%5Bsearchable%5D=true&columns%5B3%5D%5Borderable%5D=false&columns%5B3%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B3%5D%5Bsearch%5D%5Bregex%5D=false&columns%5B4%5D%5Bdata%5D=4&columns%5B4%5D%5Bname%5D=&columns%5B4%5D%5Bsearchable%5D=true&columns%5B4%5D%5Borderable%5D=false&columns%5B4%5D%5Bsearch%5D%5Bvalue%5D=&columns%5B4%5D%5Bsearch%5D%5Bregex%5D=false&order%5B0%5D%5Bcolumn%5D=0&order%5B0%5D%5Bdir%5D=asc&start=0&length=-1&search%5Bvalue%5D=&search%5Bregex%5D=false","method":"POST","mode":"cors"})
        .then(
        ''')

    def check_javascript_errors(self):
        msgs = self.driver.get_log("browser")
        for log_msg in msgs:
            print("JAVASCRIPT: %s" % log_msg['message'])
        if msgs:
            raise Exception("Javascript errors occurred")

    def get_first_gene_link(self):
        return self.driver.execute_script("""
        var ret;
        $('a').each(function(i, v){
            if($(v).attr('href') !== undefined){            
                if($(v).attr('href').slice(0,6) === '/gene/'){
                    ret = $(v).attr('href');
                    return false;
                }
            }         
        })
        return ret;
        """)

PRINT_SUBPROC = False
majiq_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

def clean_voila_files(temp_dir):
    voila_dir = os.path.join(temp_dir, 'tmp_psi')

    if os.path.exists(voila_dir):
        shutil.rmtree(voila_dir)
        os.makedirs(voila_dir)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def kill_webserver():
    try:
        check_output(('fuser', '-k', '5010/tcp'))
    except:
        print("webserver already down")

def run_majiq_build(temp_dir):

    os.environ['PYTHONPATH'] = ''
    cmd = (
           'majiq',
           'build', 'DB.gff3',
           '-c', 'settings.ini',
           '-j', '4', '-o', 'tmp_build'
           )
    # cmd = ('python3',
    #        '/home/paul/PycharmProjects/majiq/majiq/run_majiq.py',
    #        'build', 'DB.gff3',
    #        '-c', 'settings.ini',
    #        '-j', '4', '-o', 'tmp_build'
    #        )

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    for line in p.stdout:
        if PRINT_SUBPROC:
            print(line.decode().replace('\n', ''))

def run_majiq_psi(temp_dir, majiq_file, group_name):
    cmd = (
         'majiq',
         'psi', 'tmp_build/%s' % majiq_file,
         '-j', '4',
         '-o', 'tmp_psi', '-n', group_name
         )
    # cmd = ('python3',
    #        '/home/paul/PycharmProjects/majiq/majiq/run_majiq.py',
    #        'psi', 'tmp_build/%s' % majiq_file,
    #        '-j', '4',
    #        '-o', 'tmp_psi', '-n', group_name
    #        )

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    for line in p.stdout:
        if PRINT_SUBPROC:
            print(line.decode().replace('\n', ''))

def run_majiq_deltapsi(temp_dir, majiq_file1, group_name1, majiq_file2, group_name2):
    cmd = ('majiq',
           'deltapsi', '-grp1', 'tmp_build/%s' % majiq_file1,
           '-grp2', 'tmp_build/%s' % majiq_file2,
           '-j', '4',
           '-o', 'tmp_psi', '-n', group_name1, group_name2
           )
    # cmd = ('python3',
    #        '/home/paul/PycharmProjects/majiq/majiq/run_majiq.py',
    #        'deltapsi', '-grp1', 'tmp_build/%s' % majiq_file1,
    #        '-grp2', 'tmp_build/%s' % majiq_file2,
    #        '-j', '4',
    #        '-o', 'tmp_psi', '-n', group_name1, group_name2
    #        )

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    for line in p.stdout:
        if PRINT_SUBPROC:
            print(line.decode().replace('\n', ''))


def run_majiq_heterogen(temp_dir, majiq_group1_files, group_name1, majiq_group2_files, group_name2):
    cmd = ['majiq',
           'heterogen',
           '-grp1'] + ['tmp_build/%s' % f for f in majiq_group1_files] + \
           ['-grp2'] + ['tmp_build/%s' % f for f in majiq_group2_files] + \
           ['-j', '4',
           '-o', 'tmp_psi', '-n', group_name1, group_name2
            ]
    # cmd = ['python3',
    #        '/home/paul/PycharmProjects/majiq/majiq/run_majiq.py',
    #        'heterogen',
    #        '-grp1'] + ['tmp_build/%s' % f for f in majiq_group1_files] + \
    #       ['-grp2'] + ['tmp_build/%s' % f for f in majiq_group2_files] + \
    #       ['-j', '4',
    #        '-o', 'tmp_psi', '-n', group_name1, group_name2
    #        ]


    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    for line in p.stdout:
        if PRINT_SUBPROC:
            print(line.decode().replace('\n', ''))


def majiq_run(run_type, majiq_files, temp_dir):
    run_majiq_build(temp_dir)
    if run_type == 'psi':
        for file in majiq_files:
            run_majiq_psi(temp_dir, file[0], file[1])
            run_voila_tsv(temp_dir)

    elif run_type == 'deltapsi':
        run_majiq_deltapsi(temp_dir, majiq_files[0][0], majiq_files[0][1], majiq_files[1][0], majiq_files[1][1],)
        run_voila_tsv(temp_dir)
    elif run_type == 'heterogen':
        run_majiq_heterogen(temp_dir, majiq_files[0], majiq_files[1], majiq_files[2], majiq_files[3],)
        run_voila_tsv(temp_dir)

def run_voila_tsv(temp_dir):
    os.environ['PYTHONPATH'] = majiq_path
    cmd = ('python3',
           os.path.join(majiq_path, 'voila', 'run_voila.py'),
           'tsv', 'tmp_psi', 'tmp_build', '-f', '/tmp/test_debug.tsv'
           )

    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
    output = ''
    error = False
    for line in p.stdout:
        output += line.decode()
        if 'Traceback' in output:
            error = True
        if PRINT_SUBPROC:
            print(line.decode().replace('\n', ''))
    if error:
        print(output)
        assert False
    os.environ['PYTHONPATH'] = ''

def run_seq_view(outer_path, run_type, majiq_files):

    temp_dir = tempfile.mkdtemp()

    copy_tree(outer_path, temp_dir)
    #clean_voila_files(temp_dir)

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`')
    print('START TESTING %s with %s' % (run_type, str(majiq_files)))

    majiq_run(run_type, majiq_files, temp_dir)



    os.environ['PYTHONPATH'] = majiq_path

    cmd = ('python3',
           os.path.join(majiq_path, 'voila', 'run_voila.py'),
           'view', 'tmp_psi', 'tmp_build', '-p', '5010'
           )
    os.environ["PYTHONUNBUFFERED"] = "1"
    try:
        p = Popen(cmd, stdout=PIPE, stderr=STDOUT, cwd=temp_dir)
        for line in p.stdout:
            _l = line.decode().replace('\n', '')
            if PRINT_SUBPROC:
                print(_l)
            # this is not a good solution yet, for some reason we can only buffer up to the second to last line
            # so we are relying on time.sleep to wait for the web server to actually finish indexing
            if "Serving on" in _l:
                print("Voila ready")
                #time.sleep(11)
                # at this point, web server should be running
                # "creating index" is the last piece of output we seem to be able to get from subprocess
                break
        # while True:
        #     time.sleep(5)
    except KeyboardInterrupt:
        print("Interrupted")
        pass
    os.environ["PYTHONUNBUFFERED"] = ""

    return temp_dir


def test_workshop1():

    tdriver = TestWebDriver()
    temp_dir = tempfile.mkdtemp()
    try:
        majiq_files = [('workshop_Adr1.majiq', 'adr1')]
        shutil.rmtree(temp_dir)
        kill_webserver()
        temp_dir = run_seq_view(os.path.join(os.path.dirname(__file__), 'data', 'workshop_data'), 'psi', majiq_files)

        tdriver.load_page("http://localhost:5010")

        # webdriver test cases begin
        tdriver.check_javascript_errors()

        assert tdriver.get_num_records() != 0
        #assert tdriver.get_num_records() == 328

        # ------------------------


        majiq_files = [('workshop_Adr1.majiq', 'adr1'),
                       ('workshop_Adr2.majiq', 'adr2'),
                       ('workshop_Adr3.majiq', 'adr3')]

        shutil.rmtree(temp_dir)
        kill_webserver()
        temp_dir = run_seq_view(os.path.join(os.path.dirname(__file__), 'data', 'workshop_data'), 'psi', majiq_files)
        tdriver.load_page("http://localhost:5010")


        tdriver.check_javascript_errors()
        assert tdriver.get_num_records() != 0

        # ------------------------

        # majiq_files = [('workshop_Adr1.majiq', 'adr1'),
        #                ('workshop_Adr2.majiq', 'adr2'),
        #                ('workshop_Adr3.majiq', 'adr3'),
        #                ('workshop_Cer1.majiq', 'cer1'),
        #                ('workshop_Cer2.majiq', 'cer2'),
        #                ('workshop_Cer3.majiq', 'cer3')]
        # shutil.rmtree(temp_dir)
        # kill_webserver()
        # temp_dir = run_seq_view('/home/paul/PycharmProjects/majiq/test_cases/ORIGINAL_DATA', 'psi', majiq_files)
        #
        # tdriver.load_page("http://localhost:5010")
        # tdriver.check_javascript_errors()
        # assert tdriver.get_num_records() != 0

        # ------------------------

        majiq_files = [('workshop_Cer1.majiq', 'cer1'),
                       ('workshop_Adr2.majiq', 'adr2'),
                      ]
        shutil.rmtree(temp_dir)
        kill_webserver()
        temp_dir = run_seq_view(os.path.join(os.path.dirname(__file__), 'data', 'workshop_data'), 'deltapsi', majiq_files)

        tdriver.load_page("http://localhost:5010")
        tdriver.check_javascript_errors()
        assert tdriver.get_num_records() != 0
        #assert tdriver.get_num_records() == 253

        # ------------------------

        majiq_files = [('workshop_Cer1.majiq', 'workshop_Cer2.majiq', 'workshop_Cer3.majiq'),
                       'cer',
                       ('workshop_Adr1.majiq', 'workshop_Adr2.majiq', 'workshop_Adr3.majiq'),
                       'adr'
                       ]
        shutil.rmtree(temp_dir)
        kill_webserver()
        temp_dir = run_seq_view(os.path.join(os.path.dirname(__file__), 'data', 'workshop_data'), 'heterogen', majiq_files)

        tdriver.load_page("http://localhost:5010")
        tdriver.check_javascript_errors()
        assert tdriver.get_num_records() != 0
        #assert tdriver.get_num_records() == 299

        tdriver.load_page("http://localhost:5010" + tdriver.get_first_gene_link())
        tdriver.check_javascript_errors()


    except:
        print("Error running")
        print(traceback.format_exc())
    finally:
        # shut down the web server process after
        kill_webserver()
        try:
            shutil.rmtree(temp_dir)
        except:
            pass
        # remove temp files
        tdriver.end()


test_workshop1()

