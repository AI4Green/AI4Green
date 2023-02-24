import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utilities import read_yaml
from bs4 import BeautifulSoup
import requests
from urllib.request import urlopen
from shutil import copyfileobj
from datetime import datetime
import gzip
import hashlib
import ssl
context = ssl.SSLContext()

# get the compound limit from the yaml file
compound_limit = read_yaml(['compound_limit'])


def find_files(url):
    soup = BeautifulSoup(requests.get(url).text, "lxml")
    hrefs = []
    for a in soup.find_all('a'):
        hrefs.append(a['href'])
    return hrefs


def write_date_to_file():
    """Updates the date written in the download date file"""
    date_file = open('Download-Date.txt', 'w+')
    today = datetime.today()
    today = today.strftime('%Y-%m-%d')
    date_file.write(str(today))
    date_file.close()


def get_date(file):
    try:
        f = open(file, "r")
        date_last_str = f.read()
        f.close()
        date_last_obj = datetime.strptime(date_last_str, '%Y-%m-%d')
    except FileNotFoundError:
        date_last_obj = 'Unknown date'
    except ValueError:
        date_last_obj = 'Unknown date'
    return date_last_obj


def download_file(url, new_filename):
    with urlopen(url, context=context) as in_stream, \
            open(new_filename, 'wb') as out_file:
        copyfileobj(in_stream, out_file)
    print("Download Completed Successfully")


def unzip_file(filename):
    try:
        with gzip.open(filename, 'rb') as in_stream, open(filename + '_db.xml', 'wb') as out_file:
            copyfileobj(in_stream, out_file)
        print("Unzip Completed Successfully")
    except EOFError:
        print(f"The file: {filename} could not be unzipped. Delete this file and retry")
        exit()


def read_checksum_file(file):
    """Takes the checksum file and returns the checksum value as a string"""
    with open(file, "r") as f:
        text = f.read()
    all_text = text.split(" ")
    checksum = all_text[0]
    return checksum


def md5_integrity_check(file_name, checksum_file_name):
    """Checks the integrity of the download by returning the md5 for the pubchem database file
    and checks it is the same as the value given in the checksum file"""
    with open(file_name, 'rb') as file:
        data = file.read()
        md5_expected = read_checksum_file(checksum_file_name)
        md5_returned = hashlib.md5(data).hexdigest()
    if md5_expected != md5_returned:
        print("MD5 verification failed")
        exit()


def delete_temp_files():
    to_delete = ["temp_update.sqlite", "Pubchem-update", "Pubchem-update_db.xml"]
    for file in to_delete:
        try:
            os.remove(file)
        except:
            pass

