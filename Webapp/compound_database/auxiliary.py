import contextlib
import gzip
import hashlib
import os
import ssl
from datetime import datetime
from pathlib import Path
from shutil import copyfileobj
from urllib.request import urlopen

import requests
from bs4 import BeautifulSoup

context = ssl.SSLContext()

app_directory = Path(__file__).resolve().parent.parent
compound_database_dir = app_directory / "compound_database"


def find_files(url):
    soup = BeautifulSoup(requests.get(url).text, "lxml")
    hrefs = []
    for a in soup.find_all("a"):
        hrefs.append(a["href"])
    return hrefs


def write_date_to_file():
    """Updates the date written in the download date file"""
    date_file = open("Download-Date.txt", "w+")
    today = datetime.today()
    today = today.strftime("%Y-%m-%d")
    date_file.write(str(today))
    date_file.close()


def get_date(file):
    try:
        with open(file, "r") as f:
            date_last_str = f.read()
        date_last_obj = datetime.strptime(date_last_str, "%Y-%m-%d")
    except (FileNotFoundError, ValueError):
        date_last_obj = "Unknown date"
    return date_last_obj


def download_file(url, new_filename):
    with urlopen(url, context=context) as in_stream, open(
        new_filename, "wb"
    ) as out_file:
        copyfileobj(in_stream, out_file)
    print("Download Completed Successfully")


def unzip_file(infile, outfile):
    try:
        with gzip.open(infile, "rb") as f_in:
            with open(outfile, "wb") as f_out:
                copyfileobj(f_in, f_out)
    except EOFError:
        print(f"The file: {infile} could not be unzipped. Delete this file and retry")
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
    with open(file_name, "rb") as file:
        data = file.read()
        md5_expected = read_checksum_file(checksum_file_name)
        md5_returned = hashlib.md5(data).hexdigest()
    if md5_expected != md5_returned:
        print("MD5 verification failed")
        exit()


def delete_temp_files():
    to_delete = [
        "Pubchem-Database_db.xml",
        "temp_update.sqlite",
        "Pubchem-update",
        "Pubchem-update_db.xml",
    ]
    for file in to_delete:
        with contextlib.suppress(OSError):
            os.remove(compound_database_dir / file)
