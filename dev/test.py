#!/usr/bin/env python
# -*- coding: utf-8 -*-

import default_parameters
import simulation

params = default_parameters.dic

dic = {'params': params, 'data_dir': 'test', 'plotdata': True, 'verbosemode': True,
        'survey_file': 'D.txt'}

simulation.run_sim(dic)
