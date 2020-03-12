#!/usr/bin/python3


from subprocess import DEVNULL, run
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile
from shutil import unpack_archive
from threading import Thread
import logging
import platform

from novowrap.utils import accessible

log = logging.getLogger('novowrap')


def test_cmd(program, option='-v'):
    """
    Test given program and option is ok to run or not.
    Args:
        program(Path or str): program path, could be relative path if it can
        be found in $PATH or %PATH%
        option(str): option for program, usually use "-v" to show version to
        test the program
    Return:
        success(bool): success or not
    """
    test = run(f'{program} {option}', shell=True, stdout=DEVNULL,
               stderr=DEVNULL)
    success = True if test.returncode == 0 else False
    return success


def get_third_party():
    """
    Get third_party folder.
    If do not exist, create it.
    If cannot access, report.
    Return:
        success(bool): ok or not
        third_party(Path): absolute path of third_party folder
    """
    third_party = Path().home().absolute() / '.novowrap'
    success = False
    if not third_party.exists():
        log.debug(f'Create folder {third_party}')
        try:
            third_party.mkdir()
        except Exception:
            log.critical(f'Failed to create {third_party}.'
                         'Please contact the administrator.')
            return success, third_party
    if not accessible(third_party/'test', 'file'):
        log.critical(f'Failed to access {third_party}.'
                     f'Please contact the administrator.')
        return success, third_party
    success = True
    return success, third_party


def get_perl():
    """
    Linux and Mac have perl already.
    For Windows user, this function help to get perl.exe .
    Only support x86_64 or amd64 machine.
    Return:
        perl(str): perl location, empty for fail
    """
    url = ('http://strawberryperl.com/download/5.30.1.1/'
           'strawberry-perl-5.30.1.1-64bit-portable.zip')
    third_party_ok, third_party = get_third_party()
    if not third_party_ok:
        return ''
    # only Windows need it
    home_perl = third_party / 'strawberry_perl' / 'perl' / 'bin' / 'perl.exe'
    # use run instead of find_executable because the later only check if exist
    # and ignore if could run
    perl = 'perl'
    if test_cmd(perl):
        return perl
    if test_cmd(home_perl):
        return str(home_perl)
    if platform.system() != 'Windows':
        log.critical('Cannot find Perl. Please follow the link '
                     'https://www.perl.org/get.html to install.')
        return ''
    else:
        log.warning('Cannot find Perl. Try to install.')
        if '64' not in platform.machine():
            log.critical(f'Unsupport machine {platform.machine()}')
            return ''
        try:
            # file is 148mb, 148mb/3000s=~50kb/s, consider it's ok for
            # most of users
            log.info('Download may be slow. Please consider manualy install '
                     'Perl according to https://www.perl.org/get.html ')
            down = urlopen(url, timeout=3000)
        except Exception:
            log.critical('Cannot download Perl. Please try to manually '
                         ' install.')
            return ''
        zip_file = third_party / 'strawberry-perl.zip'
        with open(zip_file, 'wb') as out:
            out.write(down.read())
        folder = third_party / 'strawberry_perl'
        with ZipFile(zip_file) as z:
            z.extractall(folder)
        # fixed path in zip file, should not be wrong
        assert test_cmd(home_perl)
        return str(home_perl)


def get_novoplasty():
    """
    Ensure perl and novoplasty is available.
    Return novoplasty's path or None.
    Return:
        novoplasty(Path or None): path of perl file, None for fail
    """
    url = 'https://github.com/ndierckx/NOVOPlasty/archive/NOVOPlasty3.7.2.zip'
    third_party_ok, third_party = get_third_party()
    if not third_party_ok:
        return None
    novoplasty = (third_party / 'NOVOPlasty-NOVOPlasty3.7.2' /
                  'NOVOPlasty3.7.2.pl')
    if novoplasty.exists():
        log.debug('Found NOVOPlasty in third_party folder.')
        return novoplasty
    log.critical('Cannot find NOVOPlasty, try to download.')
    # license warning
    log.info('\tThe program assumes that user ACCEPT the license of '
             'NOVOPlasty. ')
    log.info('\tIF NOT, please terminate the running.')
    log.info('\tSee https://raw.githubusercontent.com/'
             'ndierckx/NOVOPlasty/ for details.')
    try:
        # 126kb/10s=12.6kb/s, enough for test
        _ = urlopen('https://github.com', timeout=10)
    except Exception:
        log.critical('Cannot connect to github.com.')
        log.critical('Please check your Internet connection.')
        return None
    log.info('Due to connection speed, may need minutes. Please be patient.')
    try:
        # file is ~18mb, 18mb/360s=50kb/s, consider it's ok for most of users
        down = urlopen(url, timeout=360)
    except Exception:
        log.critical('Cannot download NOVOPlasty.')
        log.critical('Please manually download it from '
                     'https://github.com/ndierckx/NOVOPlasty')
        return None
    zip_file = third_party / 'NOVOPlasty3.7.2.zip'
    with open(zip_file, 'wb') as out:
        out.write(down.read())
    # novoplasty's files have illegal character ":" which cannot be extract by
    # shutil.unpack_archive in Windows, have to use zipfile
    # windows and linux both use "/"
    with ZipFile(zip_file, 'r') as z:
        z.extractall(third_party)
    log.info(f'Got {novoplasty.stem}.')
    return novoplasty


def get_blast():
    """
    Get BLAST location.
    If BLAST was found, assume makeblastdb is found, too.
    If not found, download it.
    Return:
        ok(bool): success or not
        blast(str): blast path
    """
    third_party_ok, third_party = get_third_party()
    if not third_party_ok:
        return third_party_ok, ''
    home_blast = third_party / 'ncbi-blast-2.10.0+' / 'bin' / 'blastn'
    # in Windows, ".exe" can be omitted
    # win_home_blast = home_blast.with_name('blastn.exe')
    ok = False
    # older than 2.8.1 is buggy
    url = ('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/'
           'ncbi-blast-2.10.0+')
    urls = {'Linux': url+'-x64-linux.tar.gz',
            'Darwin': url+'-x64-macosx.tar.gz',
            'Windows': url+'-x64-win64.tar.gz'}
    blast = 'blastn'
    if test_cmd(blast, '-version'):
        ok = True
        return ok, blast
    if test_cmd(home_blast, '-version'):
        ok = True
        return ok, str(home_blast)
    log.warning('Cannot find NCBI BLAST, try to install.')
    log.info('According to Internet speed, may be slow.')
    try:
        # 50kb/10s=5kb/s, enough for test
        _ = urlopen('https://www.ncbi.nlm.nih.gov', timeout=10)
    except Exception:
        log.critical('Cannot connect to NCBI.')
        log.critical('Please check your Internet connection.')
        return ok, ''
    try:
        # file is 86-222mb, 222mb/3600s=60kb/s, consider it's ok for users all
        # over the world
        down = urlopen(urls[platform.system()], timeout=3600)
    except Exception:
        log.critical('Cannot download BLAST.')
        log.critical('Please manually download it from'
                     'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/')
        return ok, ''
    down_file = third_party / 'BLAST_2.10.0.tar.gz'
    with open(down_file, 'wb') as out:
        out.write(down.read())
    unpack_archive(down_file, third_party)
    assert test_cmd(home_blast, '-version')
    ok = True
    return ok, str(home_blast)


def third_party_main():
    """
    Use three threads to speed up.
    """
    perl = Thread(target=get_perl)
    novoplasty = Thread(target=get_novoplasty)
    blast = Thread(target=get_blast)
    perl.start()
    novoplasty.start()
    blast.start()
    perl.join()
    novoplasty.join()
    blast.join()


if __name__ == '__main__':
    third_party_main()
