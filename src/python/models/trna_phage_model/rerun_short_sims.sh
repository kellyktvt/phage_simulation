#!/bin/bash

# To kill all simulation screens:
#   screen -ls | grep sim_ | cut -d. -f1 | awk '{print $1}' | xargs -I {} screen -X -S {} quit

