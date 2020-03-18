#!/usr/bin/env python

# this file updates iconpicker.js to contain the most recent versions of icons.

import urllib.request
import yaml
import json
import re

# parse the yml file from font awesome
url = "https://raw.githubusercontent.com/FortAwesome/Font-Awesome/master/metadata/icons.yml"
f = urllib.request.urlopen(url)
myfile = f.read()
yml = yaml.load(myfile, Loader=yaml.FullLoader)

# create the icons json object, to be put into iconpicker.js
jsonf = {'icons':[]}
for d in yml: # loop through all keys of the yml dict
    for s in yml[d]['styles']: # create one object for each style (e.g. "brand", "solid") of an icon
      jsonf['icons'].append(
          {
          'title': 'fa'+s[0]+ ' fa-'+d, # create 'fab' class name for brand icons, "fas" for solid, etc.
          'searchTerms': yml[d]['search']['terms']
          }
      )
jsonf = json.dumps(jsonf, ensure_ascii=False) # prohibit ascii conversion, which breaks for special characters
f.close()

# replace icon json object inside iconpicker.js
with open('public_html/assets/js/iconpicker.js', 'r+') as f :
  jsf = f.read()
  jsf = re.sub(r'\{\n\s+icons(.*)\}',jsonf,jsf,flags=re.M|re.S)
  f.seek(0)
  f.write(jsf)
  f.truncate()
