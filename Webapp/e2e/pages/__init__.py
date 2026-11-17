"""
Page objects for the end-to-end tests.

Tests talk to these classes rather than to selectors, so when a Jinja page is
rewritten in React only its page object needs to change. Prefer role and label
locators, which survive a change of framework better than ids or class names.
"""
