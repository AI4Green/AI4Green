from selenium.common.exceptions import ElementNotInteractableException
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.select import Select
from sources.extensions import db
from utilities import basedir
from selenium import webdriver
import os
import time


def setup_selenium(self, headless='headless', maximised='maximised'):
    """Set up selenium in standard way. Requires chromedriver installation to specified path.
    Can remove headless/maximised through keyword args"""
    op = webdriver.ChromeOptions()
    if headless == 'headless':
        op.add_argument('headless')
    if maximised == 'maximised':
        op.add_argument("--start-maximized")
    chrome_dir = os.path.join(basedir, "sources/tests/test_data/chromedriver")
    self.driver = webdriver.Chrome(chrome_dir, options=op)
    self.driver.get("http://localhost:8943")


def login(self, username="PI_Test", password="PI_login"):
    """Logs in to a user account for the specified username and password - also accepts cookies"""
    try:
        self.driver.find_element_by_id('accept_cookies').click()
    except ElementNotInteractableException:
        pass
    time.sleep(1)
    clear_and_send_keys(self, 'username', username)
    clear_and_send_keys(self, 'password', password)
    self.driver.find_element_by_id('submit').click()
    time.sleep(1)
    try:
        self.driver.find_element_by_id('accept_cookies').click()
    except ElementNotInteractableException:
        pass


def select_workgroup(self, idx=1):
    select_wg = Select(self.driver.find_element_by_id("WG-select"))
    select_wg.select_by_index(idx)
    scroll_element(self, "go-to-workgroup", click=True)
    time.sleep(2)


def logout(self):
    self.driver.find_element_by_id('user-dropdown').click()
    self.driver.find_element_by_id("TopNavLoginButton").click()


def check_notifications(self):
    self.driver.find_element_by_id('user-dropdown').click()
    self.driver.find_element_by_id('notifications').click()
    time.sleep(1)


def select_workbook(self, idx=1):
    select_wb = Select(self.driver.find_element_by_id("active-workbook"))
    select_wb.select_by_index(idx)
    time.sleep(1)


def make_new_reaction(self, name='new reaction'):
    self.driver.find_element_by_id("new-reaction").click()
    time.sleep(1)
    self.driver.find_element_by_id('new-reaction-name').send_keys(name)
    self.driver.find_element_by_id('new-reaction-data-submit').click()
    time.sleep(6)


def demo_reaction(self):
    self.driver.find_element_by_id('demo-button').click()
    time.sleep(3)
    self.driver.find_element_by_id('action-button-submit').click()
    time.sleep(6)


def move_down_one_option_action_chains(self):
    actions = webdriver.ActionChains(self.driver)
    actions.send_keys(Keys.ARROW_DOWN)
    time.sleep(1)
    actions.send_keys(Keys.ENTER)
    time.sleep(1)
    actions.perform()
    time.sleep(1)


def clear_and_send_keys(self, elem, keys):
    self.driver.find_element_by_id(elem).clear()
    self.driver.find_element_by_id(elem).send_keys(keys)


def restore_db():
    dropped_tables_ls = ['Reaction', 'User', 'Person', 'WorkBook', 'WorkGroup', 'Institution', 'NovelCompound',
                         'Reaction_Reaction', 'Person_WorkBook', 'Person_WorkGroup', 'Person_WorkGroup_2',
                         'Person_WorkGroup_3', 'WorkGroup_request', 'Notification', 'WBStatusRequest',
                         'WGStatusRequest', 'NewsItem', 'ReactionNote', 'ReactionDataFile']
    for table in dropped_tables_ls:
        db.drop_table(table, with_all_data=True)
    db.disconnect()


def fill_in_physical_forms(self, phys_form_dic):
    """
    The key is the id of the compound's physical form dropdown and the value is the index of the desired
    physical form.
    Put into a looped function to try and reduce the amount of code that needs to be written
    """
    for item in phys_form_dic.items():
        print(item)
        phys_form_dropdown = self.driver.find_element_by_id(item[0])
        Select(phys_form_dropdown).select_by_index(item[1])


def add_reaction_sketcher(self, rxn_smiles, replace=False):
    """This function takes a reaction smiles string and puts it into the Marvin JS sketcher"""
    # scroll to the top of the page
    # top_page_elem = self.driver.find_element_by_id("tutorial-mode")
    # self.driver.execute_script("arguments[0].scrollIntoView(true);", top_page_elem)
    # switch to the marvin js frame
    self.driver.switch_to.frame(0)
    # click on the import button (2nd along) to open import structure popup
    elem_ls = self.driver.find_elements_by_xpath('//*[local-name() = "svg"]')
    elem_ls[1].click()
    # change the format to smiles
    format_dropdown = Select(self.driver.find_elements_by_xpath('//*[@class="gwt-ListBox"]')[0])
    format_dropdown.select_by_visible_text('SMILES')
    # find the textarea in the popup window and save the smiles string as a variable
    text_area = self.driver.find_elements_by_xpath('//*[@class="gwt-TextArea"]')[0]
    text_area.clear()
    time.sleep(0.5)
    text_area.send_keys(rxn_smiles)
    if replace:
        replace_button = self.driver.find_elements_by_xpath('//*[@class="gwt-Button mjs-marginTop"]')[1].click()
    else:
        add_button = self.driver.find_elements_by_xpath('//*[@class="gwt-Button mjs-marginTop"]')[0].click()
    self.driver.switch_to.parent_frame()
    time.sleep(1)


def sketcher_to_smiles(self):
    """ This function returns the smiles string of the reaction currently in the chemical sketcher """
    # scroll to the top of the page - prevent navbar intercepting click
    top_page_elem = self.driver.find_element_by_id("tutorial-mode")
    self.driver.execute_script("arguments[0].scrollIntoView(true);", top_page_elem)
    # switch to the marvin js frame
    self.driver.switch_to.frame(0)
    # click on the export button (3rd along) to open export structure popup
    elem_ls = self.driver.find_elements_by_xpath('//*[local-name() = "svg"]')
    elem_ls[2].click()
    # change the format to smiles
    format_dropdown = Select(self.driver.find_elements_by_xpath('//*[@class="gwt-ListBox"]')[0])
    format_dropdown.select_by_visible_text('SMILES')
    time.sleep(1)
    # find the textarea in the popup window and save the smiles string as a variable
    ids = self.driver.find_elements_by_xpath('//*[@id]')
    for idx, idd in enumerate(ids):
        if idd.tag_name == 'textarea':
            test_smiles = idd.get_attribute("value")
            return test_smiles


def scroll_element(self, elem_id, click=False):
    """Needed to avoid navbar intercepting click when scrolling back up page"""
    element = self.driver.find_element_by_id(elem_id)
    self.driver.execute_script("arguments[0].scrollIntoView(true);", element)
    if click:
        self.driver.find_element_by_id(elem_id).click()
    else:
        return element

def fill_in_reaction_table_for_demo_reaction(self):
    self.driver.find_element_by_id('js-reactant-rounded-mass1').send_keys("500")
    self.driver.find_element_by_id('js-reactant-equivalent2').send_keys("2")
    phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                     'js-product-physical-form1': 1}
    fill_in_physical_forms(self, phys_form_dic)
    self.driver.find_element_by_id('action-summary').click()
    time.sleep(8)

def fill_in_required_summary_fields(self):
    clear_and_send_keys(self, 'js-unreacted-reactant-mass', '1')
    clear_and_send_keys(self, 'js-real-product-mass', '9')
    # lock reaction
    self.driver.find_element_by_id('complete-reaction-button').click()
    time.sleep(1)


def make_and_lock_demo_reaction(self):
    make_new_reaction(self)
    demo_reaction(self)
    fill_in_reaction_table_for_demo_reaction(self)
    fill_in_required_summary_fields(self)
