#!/usr/bin/env python


"""
450k:
GSM7915081
GSM7915082
GSM7915083
GSM7915084
GSM7915282
GSM7915283
GSM7915284
GSM7915285
GSM7915288
GSM7915289

epic:
GSM7915113
GSM7915114
GSM7915157
GSM7915158
GSM7915165
GSM7915166
GSM7915167
GSM7915189
GSM7915207
GSM7915209
GSM7915234
GSM7915243
GSM7915244
"""


import urllib.request
import re
from tqdm import tqdm
import time


keys = set(['id'])

idx = {}
idx['GSM7915078'] = {'id':'glioma_Recurrence 1 [GLSS-FG-0001-R1-01M-450-AED962]'}
idx['GSM7915079'] = {'id':'glioma_Primary [GLSS-FG-0001-TP-01M-450-754721]'}
idx['GSM7915080'] = {'id':'glioma_Primary [GLSS-FG-0001-TP-02M-450-D60E26]'}
idx['GSM7915081'] = {'id':'glioma_Recurrence 1 [GLSS-FG-0002-R1-01M-450-70581A]'}
idx['GSM7915082'] = {'id':'glioma_Recurrence 2 [GLSS-FG-0002-R2-01M-450-515BA3]'}
idx['GSM7915083'] = {'id':'glioma_Primary [GLSS-FG-0002-TP-01M-450-12EF91]'}
idx['GSM7915084'] = {'id':'glioma_Primary [GLSS-FG-0002-TP-02M-450-948D90]'}
idx['GSM7915085'] = {'id':'glioma_Recurrence 1 [GLSS-FG-0003-R1-01M-450-FE0A32]'}
idx['GSM7915086'] = {'id':'glioma_Primary [GLSS-FG-0003-TP-01M-450-1B5F61]'}
idx['GSM7915087'] = {'id':'glioma_Recurrence 1 [GLSS-HF-01E6-R1-02M-EPC-EF4DDB]'}
idx['GSM7915088'] = {'id':'glioma_Primary [GLSS-HF-01E6-TP-02M-EPC-272C58]'}
idx['GSM7915089'] = {'id':'glioma_Recurrence 1 [GLSS-HF-0EAE-R1-01M-EPC-2367A2]'}
idx['GSM7915090'] = {'id':'glioma_Primary [GLSS-HF-0EAE-TP-01M-EPC-72FBD3]'}
idx['GSM7915091'] = {'id':'glioma_Recurrence 1 [GLSS-HF-1CCD-R1-01M-EPC-4BC5A1]'}
idx['GSM7915092'] = {'id':'glioma_Recurrence 2 [GLSS-HF-1CCD-R2-01M-EPC-529132]'}
idx['GSM7915093'] = {'id':'glioma_Recurrence 1 [GLSS-HF-1DDC-R1-01M-EPC-9176F0]'}
idx['GSM7915094'] = {'id':'glioma_Primary [GLSS-HF-1DDC-TP-01M-EPC-BC0FE1]'}
idx['GSM7915095'] = {'id':'glioma_Recurrence 1 [GLSS-HF-1E83-R1-01M-EPC-362E33]'}
idx['GSM7915096'] = {'id':'glioma_Recurrence 2 [GLSS-HF-1E83-R2-01M-EPC-2B7636]'}
idx['GSM7915097'] = {'id':'glioma_Primary [GLSS-HF-1E83-TP-01M-EPC-531D30]'}
idx['GSM7915098'] = {'id':'glioma_Recurrence 1 [GLSS-HF-2548-R1-02M-EPC-280C13]'}
idx['GSM7915099'] = {'id':'glioma_Primary [GLSS-HF-2548-TP-02M-EPC-F65D55]'}
idx['GSM7915100'] = {'id':'glioma_Recurrence 1 [GLSS-HF-2DF0-R1-01M-EPC-82F06C]'}
idx['GSM7915101'] = {'id':'glioma_Primary [GLSS-HF-2DF0-TP-01M-EPC-3CEC1B]'}
idx['GSM7915102'] = {'id':'glioma_Recurrence 1 [GLSS-HF-3050-R1-01M-EPC-914E6C]'}
idx['GSM7915103'] = {'id':'glioma_Recurrence 2 [GLSS-HF-3050-R2-01M-EPC-29ACCE]'}
idx['GSM7915104'] = {'id':'glioma_Primary [GLSS-HF-3050-TP-01M-EPC-ACD195]'}
idx['GSM7915105'] = {'id':'glioma_Recurrence 1 [GLSS-HF-3118-R1-01M-EPC-7645C3]'}
idx['GSM7915106'] = {'id':'glioma_Primary [GLSS-HF-3118-TP-01M-EPC-F7E790]'}
idx['GSM7915107'] = {'id':'glioma_Recurrence 1 [GLSS-HF-3162-R1-02M-EPC-109017]'}
idx['GSM7915108'] = {'id':'glioma_Primary [GLSS-HF-3162-TP-02M-EPC-4A333F]'}
idx['GSM7915109'] = {'id':'glioma_Recurrence 1 [GLSS-HF-38DA-R1-01M-EPC-C52925]'}
idx['GSM7915110'] = {'id':'glioma_Primary [GLSS-HF-38DA-TP-01M-EPC-8199DB]'}
idx['GSM7915111'] = {'id':'glioma_Recurrence 1 [GLSS-HF-409D-R1-01M-EPC-2AACCA]'}
idx['GSM7915112'] = {'id':'glioma_Primary [GLSS-HF-409D-TP-01M-EPC-62ABD1]'}
idx['GSM7915113'] = {'id':'glioma_Recurrence 1 [GLSS-HF-4B46-R1-01M-EPC-455B6E]'}
idx['GSM7915114'] = {'id':'glioma_Primary [GLSS-HF-4B46-TP-01M-EPC-A920DE]'}
idx['GSM7915115'] = {'id':'glioma_Recurrence 1 [GLSS-HF-4D9A-R1-01M-EPC-AF1FB6]'}
idx['GSM7915116'] = {'id':'glioma_Primary [GLSS-HF-4D9A-TP-01M-EPC-C80747]'}
idx['GSM7915117'] = {'id':'glioma_Recurrence 1 [GLSS-HF-4F0A-R1-01M-EPC-6642F6]'}
idx['GSM7915118'] = {'id':'glioma_Recurrence 2 [GLSS-HF-4F0A-R2-01M-EPC-E4FA81]'}
idx['GSM7915119'] = {'id':'glioma_Primary [GLSS-HF-4F0A-TP-01M-EPC-840C2F]'}
idx['GSM7915120'] = {'id':'glioma_Primary [GLSS-HF-4FBF-R1-01M-EPC-C3406D]'}
idx['GSM7915121'] = {'id':'glioma_Recurrence 1 [GLSS-HF-4FBF-R2-01M-EPC-1CEBCD]'}
idx['GSM7915122'] = {'id':'glioma_Recurrence 2 [GLSS-HF-50F3-R2-01M-EPC-6753AA]'}
idx['GSM7915123'] = {'id':'glioma_Primary [GLSS-HF-50F3-TP-01M-EPC-93C210]'}
idx['GSM7915124'] = {'id':'glioma_Recurrence 2 [GLSS-HF-53B9-R2-01M-EPC-A0A5F3]'}
idx['GSM7915125'] = {'id':'glioma_Primary [GLSS-HF-53B9-TP-01M-EPC-2ED903]'}
idx['GSM7915126'] = {'id':'glioma_Recurrence 1 [GLSS-HF-57AE-R1-01M-EPC-5AC507]'}
idx['GSM7915127'] = {'id':'glioma_Recurrence 2 [GLSS-HF-57AE-R2-01M-EPC-8ADA14]'}
idx['GSM7915128'] = {'id':'glioma_Primary [GLSS-HF-57AE-TP-01M-EPC-3A1C51]'}
idx['GSM7915129'] = {'id':'glioma_Recurrence 1 [GLSS-HF-630A-R1-01M-EPC-7532E2]'}
idx['GSM7915130'] = {'id':'glioma_Primary [GLSS-HF-630A-TP-01M-EPC-92AEC9]'}
idx['GSM7915131'] = {'id':'glioma_Recurrence 1 [GLSS-HF-6504-R1-01M-EPC-252EE9]'}
idx['GSM7915132'] = {'id':'glioma_Primary [GLSS-HF-6504-TP-01M-EPC-E06884]'}
idx['GSM7915133'] = {'id':'glioma_Recurrence 1 [GLSS-HF-6658-R1-01M-EPC-0EDB01]'}
idx['GSM7915134'] = {'id':'glioma_Primary [GLSS-HF-6658-TP-01M-EPC-EFC2A0]'}
idx['GSM7915135'] = {'id':'glioma_Recurrence 2 [GLSS-HF-753F-R2-01M-EPC-2AF3BB]'}
idx['GSM7915136'] = {'id':'glioma_Recurrence 3 [GLSS-HF-753F-R3-01M-EPC-C4C6B3]'}
idx['GSM7915137'] = {'id':'glioma_Recurrence 4 [GLSS-HF-753F-R4-01M-EPC-76AAC7]'}
idx['GSM7915138'] = {'id':'glioma_Recurrence 1 [GLSS-HF-7CAC-R1-01M-EPC-27477B]'}
idx['GSM7915139'] = {'id':'glioma_Primary [GLSS-HF-7CAC-TP-01M-EPC-B10748]'}
idx['GSM7915140'] = {'id':'glioma_Recurrence 1 [GLSS-HF-8438-R1-01M-EPC-885AEF]'}
idx['GSM7915141'] = {'id':'glioma_Recurrence 2 [GLSS-HF-8438-R2-01M-EPC-A9D6C2]'}
idx['GSM7915142'] = {'id':'glioma_Primary [GLSS-HF-8438-TP-01M-EPC-37E08B]'}
idx['GSM7915143'] = {'id':'glioma_Recurrence 1 [GLSS-HF-887A-R1-01M-EPC-9BA269]'}
idx['GSM7915144'] = {'id':'glioma_Primary [GLSS-HF-887A-TP-01M-EPC-89FB60]'}
idx['GSM7915145'] = {'id':'glioma_Recurrence 1 [GLSS-HF-891E-R1-01M-EPC-DD5035]'}
idx['GSM7915146'] = {'id':'glioma_Recurrence 2 [GLSS-HF-891E-R2-01M-EPC-9FE93D]'}
idx['GSM7915147'] = {'id':'glioma_Primary [GLSS-HF-891E-TP-01M-EPC-4FE57C]'}
idx['GSM7915148'] = {'id':'glioma_Recurrence 1 [GLSS-HF-8FCD-R1-01M-EPC-7F42B9]'}
idx['GSM7915149'] = {'id':'glioma_Primary [GLSS-HF-8FCD-TP-01M-EPC-A16FA0]'}
idx['GSM7915150'] = {'id':'glioma_Recurrence 1 [GLSS-HF-9A7A-R1-01M-EPC-7D99F6]'}
idx['GSM7915151'] = {'id':'glioma_Recurrence 2 [GLSS-HF-9A7A-R2-01M-EPC-A077DA]'}
idx['GSM7915152'] = {'id':'glioma_Primary [GLSS-HF-9A7A-TP-01M-EPC-E7F59A]'}
idx['GSM7915153'] = {'id':'glioma_Recurrence 1 [GLSS-HF-AD84-R1-01M-EPC-FAD820]'}
idx['GSM7915154'] = {'id':'glioma_Primary [GLSS-HF-AD84-TP-01M-EPC-34B931]'}
idx['GSM7915155'] = {'id':'glioma_Recurrence 2 [GLSS-HF-B30B-R2-01M-EPC-75D063]'}
idx['GSM7915156'] = {'id':'glioma_Primary [GLSS-HF-B30B-TP-01M-EPC-077905]'}
idx['GSM7915157'] = {'id':'glioma_Recurrence 1 [GLSS-HF-B3DE-R1-01M-EPC-DDF7EA]'}
idx['GSM7915158'] = {'id':'glioma_Primary [GLSS-HF-B3DE-TP-01M-EPC-F89A8B]'}
idx['GSM7915159'] = {'id':'glioma_Recurrence 2 [GLSS-HF-B92C-R2-01M-EPC-58A139]'}
idx['GSM7915160'] = {'id':'glioma_Primary [GLSS-HF-B92C-TP-01M-EPC-C2F257]'}
idx['GSM7915161'] = {'id':'glioma_Recurrence 1 [GLSS-HF-B972-R1-01M-EPC-9AE184]'}
idx['GSM7915162'] = {'id':'glioma_Primary [GLSS-HF-B972-TP-01M-EPC-5A550B]'}
idx['GSM7915163'] = {'id':'glioma_Recurrence 3 [GLSS-HF-BCDE-R3-01M-EPC-ECF8B2]'}
idx['GSM7915164'] = {'id':'glioma_Primary [GLSS-HF-BCDE-TP-01M-EPC-064E29]'}
idx['GSM7915165'] = {'id':'glioma_Recurrence 1 [GLSS-HF-CB66-R1-01M-EPC-D91328]'}
idx['GSM7915166'] = {'id':'glioma_Recurrence 2 [GLSS-HF-CB66-R2-01M-EPC-1522F8]'}
idx['GSM7915167'] = {'id':'glioma_Primary [GLSS-HF-CB66-TP-01M-EPC-256484]'}
idx['GSM7915168'] = {'id':'glioma_Recurrence 1 [GLSS-HF-D48E-R1-01M-EPC-3D15C5]'}
idx['GSM7915169'] = {'id':'glioma_Primary [GLSS-HF-D48E-TP-01M-EPC-9D2248]'}
idx['GSM7915170'] = {'id':'glioma_Recurrence 1 [GLSS-HF-DCED-R1-01M-EPC-8151C0]'}
idx['GSM7915171'] = {'id':'glioma_Recurrence 2 [GLSS-HF-DCED-R2-01M-EPC-114B69]'}
idx['GSM7915172'] = {'id':'glioma_Recurrence 2 [GLSS-HF-DE05-R1-01M-EPC-ADC829]'}
idx['GSM7915173'] = {'id':'glioma_Primary [GLSS-HF-DE05-TP-01M-EPC-7B687D]'}
idx['GSM7915174'] = {'id':'glioma_Recurrence 1 [GLSS-HF-DF35-R1-01M-EPC-359F1E]'}
idx['GSM7915175'] = {'id':'glioma_Primary [GLSS-HF-DF35-TP-01M-EPC-5AA41D]'}
idx['GSM7915176'] = {'id':'glioma_Recurrence 1 [GLSS-HF-E362-R1-01M-EPC-4BEF87]'}
idx['GSM7915177'] = {'id':'glioma_Primary [GLSS-HF-E362-TP-01M-EPC-F32415]'}
idx['GSM7915178'] = {'id':'glioma_Recurrence 1 [GLSS-HF-E6BA-R1-01M-EPC-10A328]'}
idx['GSM7915179'] = {'id':'glioma_Recurrence 2 [GLSS-HF-E6BA-R2-01M-EPC-F41EBF]'}
idx['GSM7915180'] = {'id':'glioma_Recurrence 1 [GLSS-HF-EDB1-R1-01M-EPC-5D5CEF]'}
idx['GSM7915181'] = {'id':'glioma_Primary [GLSS-HF-EDB1-TP-01M-EPC-204754]'}
idx['GSM7915182'] = {'id':'glioma_Recurrence 1 [GLSS-HF-EE74-R1-01M-EPC-5AD560]'}
idx['GSM7915183'] = {'id':'glioma_Primary [GLSS-HF-EE74-TP-01M-EPC-891CEB]'}
idx['GSM7915184'] = {'id':'glioma_Recurrence 1 [GLSS-HF-EE77-R1-01M-EPC-08C083]'}
idx['GSM7915185'] = {'id':'glioma_Recurrence 2 [GLSS-HF-EE77-R2-01M-EPC-D9FE3D]'}
idx['GSM7915186'] = {'id':'glioma_Recurrence 3 [GLSS-HF-EE77-R3-01M-EPC-6AF605]'}
idx['GSM7915187'] = {'id':'glioma_Recurrence 1 [GLSS-HF-F922-R1-01M-EPC-00A2E0]'}
idx['GSM7915188'] = {'id':'glioma_Recurrence 2 [GLSS-HF-F922-R2-01M-EPC-D67141]'}
idx['GSM7915189'] = {'id':'glioma_Primary [GLSS-HF-F922-TP-01M-EPC-BED05A]'}
idx['GSM7915190'] = {'id':'glioma_Recurrence 1 [GLSS-HK-0001-R1-01M-EPC-CBF116]'}
idx['GSM7915191'] = {'id':'glioma_Primary [GLSS-HK-0001-TP-01M-EPC-966D7D]'}
idx['GSM7915192'] = {'id':'glioma_Recurrence 1 [GLSS-HK-0002-R1-01M-EPC-DBE704]'}
idx['GSM7915193'] = {'id':'glioma_Primary [GLSS-HK-0002-TP-01M-EPC-139865]'}
idx['GSM7915194'] = {'id':'glioma_Recurrence 1 [GLSS-HK-0003-R1-01M-EPC-AC36CA]'}
idx['GSM7915195'] = {'id':'glioma_Primary [GLSS-HK-0003-TP-01M-EPC-6834EA]'}
idx['GSM7915196'] = {'id':'glioma_Recurrence 1 [GLSS-LU-00B9-R1-01M-450-2D8DD7]'}
idx['GSM7915197'] = {'id':'glioma_Primary [GLSS-LU-00B9-TP-01M-450-8488BF]'}
idx['GSM7915198'] = {'id':'glioma_Recurrence 1 [GLSS-LU-00C4-R1-01M-450-F506FB]'}
idx['GSM7915199'] = {'id':'glioma_Primary [GLSS-LU-00C4-TP-01M-450-73E71B]'}
idx['GSM7915200'] = {'id':'glioma_Recurrence 1 [GLSS-LU-00C7-R1-01M-450-AC08D4]'}
idx['GSM7915201'] = {'id':'glioma_Primary [GLSS-LU-00C7-TP-01M-450-9EDFD9]'}
idx['GSM7915202'] = {'id':'glioma_Recurrence 1 [GLSS-LU-0B13-R1-01M-450-2D305E]'}
idx['GSM7915203'] = {'id':'glioma_Primary [GLSS-LU-0B13-TP-01M-450-F3AA51]'}
idx['GSM7915204'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0003-R1-01M-EPC-D7E196]'}
idx['GSM7915205'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0003-R2-01M-EPC-D5E346]'}
idx['GSM7915206'] = {'id':'glioma_Primary [GLSS-LX-0003-TP-01M-EPC-13A6E5]'}
idx['GSM7915207'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0015-R1-01M-EPC-25434C]'}
idx['GSM7915208'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0015-R2-01M-EPC-719AA5]'}
idx['GSM7915209'] = {'id':'glioma_Primary [GLSS-LX-0015-TP-01M-EPC-5C276E]'}
idx['GSM7915210'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0019-R1-01M-EPC-C5E3F5]'}
idx['GSM7915211'] = {'id':'glioma_Primary [GLSS-LX-0019-TP-01M-EPC-E33533]'}
idx['GSM7915212'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0040-R1-01M-EPC-F67FA3]'}
idx['GSM7915213'] = {'id':'glioma_Primary [GLSS-LX-0040-TP-01M-EPC-C3889A]'}
idx['GSM7915214'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0070-R1-01M-EPC-D362E1]'}
idx['GSM7915215'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0070-R2-01M-EPC-1DA5EB]'}
idx['GSM7915216'] = {'id':'glioma_Primary [GLSS-LX-0070-TP-01M-EPC-E40053]'}
idx['GSM7915217'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0072-R1-01M-EPC-2E53DF]'}
idx['GSM7915218'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0072-R2-01M-EPC-66DAB8]'}
idx['GSM7915219'] = {'id':'glioma_Primary [GLSS-LX-0072-TP-01M-EPC-57D2CA]'}
idx['GSM7915220'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0083-R1-01M-EPC-504A59]'}
idx['GSM7915221'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0083-R2-01M-EPC-07FFAE]'}
idx['GSM7915222'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0114-R1-01M-EPC-06697C]'}
idx['GSM7915223'] = {'id':'glioma_Primary [GLSS-LX-0114-TP-01M-EPC-7B03DA]'}
idx['GSM7915224'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0123-R1-01M-EPC-76C250]'}
idx['GSM7915225'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0123-R2-01M-EPC-CCC2A5]'}
idx['GSM7915226'] = {'id':'glioma_Primary [GLSS-LX-0123-TP-01M-EPC-B9E8D3]'}
idx['GSM7915227'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0132-R1-01M-EPC-728873]'}
idx['GSM7915228'] = {'id':'glioma_Primary [GLSS-LX-0132-TP-01M-EPC-045261]'}
idx['GSM7915229'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0134-R1-01M-EPC-2D1682]'}
idx['GSM7915230'] = {'id':'glioma_Primary [GLSS-LX-0134-TP-01M-EPC-CDBA63]'}
idx['GSM7915231'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0164-R1-01M-EPC-C604DC]'}
idx['GSM7915232'] = {'id':'glioma_Primary [GLSS-LX-0164-TP-01M-EPC-D09F60]'}
idx['GSM7915233'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0170-R1-01M-EPC-9558D1]'}
idx['GSM7915234'] = {'id':'glioma_Primary [GLSS-LX-0170-TP-01M-EPC-8CA2CF]'}
idx['GSM7915235'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0176-R1-01M-EPC-910551]'}
idx['GSM7915236'] = {'id':'glioma_Primary [GLSS-LX-0176-TP-01M-EPC-918026]'}
idx['GSM7915237'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0178-R1-01M-EPC-81AFEA]'}
idx['GSM7915238'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0178-R2-01M-EPC-8AAC50]'}
idx['GSM7915239'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0200-R1-01M-EPC-BB384B]'}
idx['GSM7915240'] = {'id':'glioma_Primary [GLSS-LX-0200-TP-01M-EPC-6AF76B]'}
idx['GSM7915241'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0215-R1-01M-EPC-EA4FD1]'}
idx['GSM7915242'] = {'id':'glioma_Primary [GLSS-LX-0215-TP-01M-EPC-A8603E]'}
idx['GSM7915243'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0242-R1-01M-EPC-D53E3A]'}
idx['GSM7915244'] = {'id':'glioma_Primary [GLSS-LX-0242-TP-01M-EPC-EB0091]'}
idx['GSM7915245'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0267-R1-01M-EPC-76260E]'}
idx['GSM7915246'] = {'id':'glioma_Primary [GLSS-LX-0267-TP-01M-EPC-61490B]'}
idx['GSM7915247'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0299-R1-01M-EPC-B573EC]'}
idx['GSM7915248'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0299-R2-01M-EPC-C17044]'}
idx['GSM7915249'] = {'id':'glioma_Primary [GLSS-LX-0299-TP-01M-EPC-1F29BA]'}
idx['GSM7915250'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0304-R1-01M-EPC-545F1F]'}
idx['GSM7915251'] = {'id':'glioma_Primary [GLSS-LX-0304-TP-01M-EPC-BFEB6C]'}
idx['GSM7915252'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0344-R1-01M-EPC-ECAC06]'}
idx['GSM7915253'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0344-R2-01M-EPC-3C63F0]'}
idx['GSM7915254'] = {'id':'glioma_Primary [GLSS-LX-0344-TP-01M-EPC-85C722]'}
idx['GSM7915255'] = {'id':'glioma_Recurrence 1 [GLSS-LX-0357-R1-01M-EPC-1BF185]'}
idx['GSM7915256'] = {'id':'glioma_Recurrence 2 [GLSS-LX-0357-R2-01M-EPC-6F8AB0]'}
idx['GSM7915257'] = {'id':'glioma_Recurrence 3 [GLSS-LX-0357-R3-01M-EPC-A63D4B]'}
idx['GSM7915258'] = {'id':'glioma_Recurrence 1 [TCGA-06-0125-R1-11M-450-7VHZXA]'}
idx['GSM7915259'] = {'id':'glioma_Primary [TCGA-06-0125-TP-01M-450-XC6ZJF]'}
idx['GSM7915260'] = {'id':'glioma_Recurrence 1 [TCGA-06-0152-R1-01M-450-OWMIQB]'}
idx['GSM7915261'] = {'id':'glioma_Primary [TCGA-06-0152-TP-02M-450-GGWYUC]'}
idx['GSM7915262'] = {'id':'glioma_Recurrence 1 [TCGA-06-0171-R1-11M-450-KRJEHO]'}
idx['GSM7915263'] = {'id':'glioma_Primary [TCGA-06-0171-TP-02M-450-LMVVPP]'}
idx['GSM7915264'] = {'id':'glioma_Recurrence 1 [TCGA-06-0190-R1-01M-450-BR0QI9]'}
idx['GSM7915265'] = {'id':'glioma_Primary [TCGA-06-0190-TP-01M-450-YDSMF9]'}
idx['GSM7915266'] = {'id':'glioma_Recurrence 1 [TCGA-06-0210-R1-01M-450-DNWLII]'}
idx['GSM7915267'] = {'id':'glioma_Primary [TCGA-06-0210-TP-01M-450-IST5AV]'}
idx['GSM7915268'] = {'id':'glioma_Recurrence 1 [TCGA-06-0211-R1-02M-450-PML9WT]'}
idx['GSM7915269'] = {'id':'glioma_Primary [TCGA-06-0211-TP-01M-450-7QXLC9]'}
idx['GSM7915270'] = {'id':'glioma_Recurrence 1 [TCGA-06-0221-R1-11M-450-BV33CH]'}
idx['GSM7915271'] = {'id':'glioma_Primary [TCGA-06-0221-TP-01M-450-Z6T93Z]'}
idx['GSM7915272'] = {'id':'glioma_Recurrence 1 [TCGA-14-0736-R1-01M-450-KIHOXS]'}
idx['GSM7915273'] = {'id':'glioma_Primary [TCGA-14-0736-TP-01M-450-1VSUVF]'}
idx['GSM7915274'] = {'id':'glioma_Recurrence 1 [TCGA-14-1402-R1-01M-450-MKAY9C]'}
idx['GSM7915275'] = {'id':'glioma_Primary [TCGA-14-1402-TP-01M-450-Z1WJEV]'}
idx['GSM7915276'] = {'id':'glioma_Recurrence 1 [TCGA-19-0957-R1-11M-450-IWFIBJ]'}
idx['GSM7915277'] = {'id':'glioma_Primary [TCGA-19-0957-TP-01M-450-58TSX9]'}
idx['GSM7915278'] = {'id':'glioma_Recurrence 1 [TCGA-19-1389-R1-21M-450-9WZ33C]'}
idx['GSM7915279'] = {'id':'glioma_Primary [TCGA-19-1389-TP-01M-450-CSRJQZ]'}
idx['GSM7915280'] = {'id':'glioma_Recurrence 1 [TCGA-19-4065-R1-11M-450-T2ZZWK]'}
idx['GSM7915281'] = {'id':'glioma_Primary [TCGA-19-4065-TP-01M-450-7DXWB2]'}
idx['GSM7915282'] = {'id':'glioma_Recurrence 1 [TCGA-DH-A669-R1-11M-450-W7XM5M]'}
idx['GSM7915283'] = {'id':'glioma_Primary [TCGA-DH-A669-TP-12M-450-QDG6WF]'}
idx['GSM7915284'] = {'id':'glioma_Recurrence 1 [TCGA-DU-5870-R1-12M-450-T5T4I6]'}
idx['GSM7915285'] = {'id':'glioma_Primary [TCGA-DU-5870-TP-11M-450-FOTKBT]'}
idx['GSM7915286'] = {'id':'glioma_Recurrence 1 [TCGA-DU-5872-R1-21M-450-LXD7NY]'}
idx['GSM7915287'] = {'id':'glioma_Primary [TCGA-DU-5872-TP-11M-450-EGJX4P]'}
idx['GSM7915288'] = {'id':'glioma_Recurrence 1 [TCGA-DU-6397-R1-12M-450-U8QZUQ]'}
idx['GSM7915289'] = {'id':'glioma_Primary [TCGA-DU-6397-TP-11M-450-5LIWE1]'}
idx['GSM7915290'] = {'id':'glioma_Recurrence 1 [TCGA-DU-6404-R1-21M-450-OEKA4K]'}
idx['GSM7915291'] = {'id':'glioma_Recurrence 2 [TCGA-DU-6404-R2-11M-450-8SL9O5]'}
idx['GSM7915292'] = {'id':'glioma_Primary [TCGA-DU-6404-TP-11M-450-1WUO60]'}
idx['GSM7915293'] = {'id':'glioma_Recurrence 1 [TCGA-DU-7304-R1-12M-450-ZDRTKO]'}
idx['GSM7915294'] = {'id':'glioma_Primary [TCGA-DU-7304-TP-12M-450-VQITWQ]'}
idx['GSM7915295'] = {'id':'glioma_Recurrence 1 [TCGA-FG-5963-R1-12M-450-KICV7S]'}
idx['GSM7915296'] = {'id':'glioma_Primary [TCGA-FG-5963-TP-11M-450-67KBIY]'}
idx['GSM7915297'] = {'id':'glioma_Recurrence 1 [TCGA-FG-5965-R1-11M-450-T3VSRV]'}
idx['GSM7915298'] = {'id':'glioma_Recurrence 2 [TCGA-FG-5965-R2-11M-450-YDZD5D]'}
idx['GSM7915299'] = {'id':'glioma_Primary [TCGA-FG-5965-TP-11M-450-4K77U2]'}
idx['GSM7915300'] = {'id':'glioma_Recurrence 1 [TCGA-FG-6691-R1-01M-450-211EB3]'}
idx['GSM7915301'] = {'id':'glioma_Primary [TCGA-FG-6691-TP-01M-450-CA46A7]'}
idx['GSM7915302'] = {'id':'glioma_Recurrence 1 [TCGA-FG-A4MT-R1-11M-450-P2ILGT]'}
idx['GSM7915303'] = {'id':'glioma_Primary [TCGA-FG-A4MT-TP-11M-450-6RSSPS]'}
idx['GSM7915304'] = {'id':'glioma_Recurrence 1 [TCGA-TM-A7CF-R1-11M-450-IVILY5]'}
idx['GSM7915305'] = {'id':'glioma_Primary [TCGA-TM-A7CF-TP-11M-450-A7B9ED]'}
idx['GSM7915306'] = {'id':'glioma_Recurrence 1 [TCGA-TQ-A7RK-R1-11M-450-0L00IF]'}
idx['GSM7915307'] = {'id':'glioma_Recurrence 2 [TCGA-TQ-A7RK-R2-11M-450-8ZHO2D]'}
idx['GSM7915308'] = {'id':'glioma_Primary [TCGA-TQ-A7RK-TP-11M-450-61T2CA]'}
idx['GSM7915309'] = {'id':'glioma_Recurrence 1 [TCGA-TQ-A7RR-R1-01M-450-01F584]'}
idx['GSM7915310'] = {'id':'glioma_Primary [TCGA-TQ-A7RR-TP-01M-450-6C850C]'}
idx['GSM7915311'] = {'id':'glioma_Recurrence 1 [TCGA-TQ-A7RV-R1-11M-450-3ETTB6]'}
idx['GSM7915312'] = {'id':'glioma_Primary [TCGA-TQ-A7RV-TP-21M-450-TV0QJ2]'}
idx['GSM7915313'] = {'id':'glioma_Recurrence 1 [TCGA-TQ-A7RW-R1-01M-450-66E0A6]'}
idx['GSM7915314'] = {'id':'glioma_Primary [TCGA-TQ-A7RW-TP-01M-450-0F3063]'}
idx['GSM7915315'] = {'id':'glioma_Recurrence 1 [TCGA-TQ-A8XE-R1-11M-450-WLFHXP]'}
idx['GSM7915316'] = {'id':'glioma_Primary [TCGA-TQ-A8XE-TP-11M-450-OD11JS]'}
idx['GSM7915317'] = {'id':'glioma_Recurrence 1 [TCGA-DU-6407-R1-02M-EPC-A7F9D5]'}
idx['GSM7915318'] = {'id':'glioma_Recurrence 2 [TCGA-DU-6407-R2-02M-EPC-F64EF7]'}
idx['GSM7915319'] = {'id':'glioma_Primary [TCGA-DU-6407-TP-02M-EPC-4E5F5E]'}
idx['GSM7915320'] = {'id':'glioma_Recurrence 1 [TCGA-DU-6407-R1-12M-450-V1ZOL5]'}
idx['GSM7915321'] = {'id':'glioma_Recurrence 2 [TCGA-DU-6407-R2-11M-450-DU6R11]'}
idx['GSM7915322'] = {'id':'glioma_Primary [TCGA-DU-6407-TP-13M-450-2G6KI4]'}




for geo_id in tqdm(idx.keys()):
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + geo_id
    print(url)
    
    response = urllib.request.urlopen(url)
    data = response.read()      # a `bytes` object
    text = data.decode('utf-8') # a `str`; this step can't be used if data is binary
    
    # Characteristics
    parse = text.split("Characteristics")[1]
    parse = parse.split("Treatment pro")[0]
    parse = parse.replace("<br>", "\n")
    parse = parse.replace('<td style="text-align: justify">',"")
    parse = parse.replace('<tr valign="top">',"")
    parse = parse.replace('<td nowrap>',"")
    parse = parse.replace("</tr>","")
    parse = parse.replace("</td>","")
    parse = parse.strip()
    parse = parse.split("\n")
    parse = {_.split(": ")[0]: _.split(": ")[1] for _ in parse}
    
    for key in parse.keys():
      keys.add(key)
      idx[geo_id][key] = parse[key]
    
    time.sleep(1.5)


#print(keys)
#print("\n\n---\n\n")
#print(idx)


with open("output/tables/GSM7915078.metadata.txt", "w") as fh:
  fh.write("sid")
  for key in keys:
    fh.write("\t" + key)
  fh.write("\n")
  
  for sid in idx.keys():
    fh.write(sid)
    
    for key in keys:
      if key in idx[sid]:
        fh.write("\t"+idx[sid][key])
      else:
        print("-")
    
    fh.write("\n")

