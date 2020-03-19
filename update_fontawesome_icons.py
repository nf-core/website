#!/usr/bin/env python

# Fetches the latest list of the FontAwesome icons and builds a
# JSON file with identifiers and search terms.

import urllib.request
import json

# Note that the version of FontAwesome must match the one imported in header.php
fontawesome_version = '5.12.1'

# parse the JSON file from Font Awesome
url = "https://raw.githubusercontent.com/FortAwesome/Font-Awesome/{}/metadata/icons.json".format(fontawesome_version)
with urllib.request.urlopen(url) as f:
    fa_icons = json.loads(f.read().decode())

# create the icons json object, to be put into iconpicker.js
icons = {};
for d in fa_icons: # loop through all keys of the fa_icons dict
    for s in fa_icons[d]['styles']: # create one object for each style (e.g. "brand", "solid") of an icon
      icon_class = 'fa'+s[0]+ ' fa-'+d # create 'fab' class name for brand icons, "fas" for solid, etc.
      search_terms = ' '.join(fa_icons[d]['search']['terms'])
      icons[icon_class] = search_terms

# Write to a JSON file for the website
with open('public_html/assets/js/fa-icons.json', 'w') as fh :
    json.dump(icons, fh, ensure_ascii=False) # prohibit ascii conversion, which breaks for special characters
