import os, sys
import numpy as np
from astropy.table import Table

from selenium import webdriver
from selenium.webdriver.support.select import Select


def WISE_download_command(input_file, output_file, catalog='catwise_2020', radius=2):
    command = 'curl -F filename=@%s'%(input_file)\
            + ' -F catalog=%s'%(catalog)\
            + ' -F spatial=upload'\
            + ' -F uradius=%.2f'%(radius)\
            + ' -F outfmt=1'\
            + ' -F selcols="source_name,source_id,ra,dec,w1mpro,w1sigmpro,w2mpro,w2sigmpro"'\
            + ' "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"'\
            + ' -o %s'%(output_file)\

    return command

def script_CatWISE(tbl, ra_col, dec_col, script_name, outdir,\
                   n_onematch=200000, overwrite=False):
    script = open(script_name, 'w')

    Niter = int(len(tbl)/n_onematch)

    for niter in range(Niter+1):

        input_file = os.path.abspath(outdir+'/upload/coord%d.tbl'%(niter))
        output_file = os.path.abspath(outdir+'/download/match%d.tbl'%(niter))

        if (not overwrite) and os.path.exists(output_file):
            continue

        start = niter * n_onematch
        end = np.min([(niter+1)*n_onematch, len(tbl)])
        print('Generating coordinate list for objects %d to %d'%(start, end-1))
        coord_tbl = tbl[ra_col, dec_col][start:end]
        coord_tbl.rename_column(ra_col, 'ra')
        coord_tbl.rename_column(dec_col, 'dec')
        coord_tbl.write(outdir+'/upload/coord%d.tbl'%(niter),\
                        format='ascii.ipac', overwrite=True)

        cmd = WISE_download_command(os.path.basename(input_file), output_file)
        script.write(cmd)
        script.write('\n')

    script.close()

def combine_result(ntask, output):
    alltbl = []
    for root, namedir, names in os.walk('./combine'):
        for name in names:
            if '%d_matched.tbl'%(ntask) in name:
                thistbl = Table.read(root + '/' + name, format='ascii.ipac')
                alltbl.append(thistbl)
    alltbl = vstack(alltbl)
    alltbl.write(output)

def download_CatWISE(tbl, ra_col, dec_col, outdir,\
                    script_name='./download_script.txt',
                    overwrite=False):
#    if not os.path.exists(outdir):
    outdir = os.path.abspath(outdir)

    os.system('mkdir -p %s'%outdir)
    os.system('mkdir -p %s/upload'%outdir)
    os.system('mkdir -p %s/download'%outdir)

    cwd = os.getcwd()
    os.chdir('%s/upload'%outdir)
    script_CatWISE(tbl, ra_col, dec_col, script_name, outdir,\
                  overwrite=overwrite)
    os.system('chmod a+x %s'%script_name)
    os.system(script_name)

    os.chdir(cwd)

def VHS_CrossID(filename, prog='110'):
    driver = webdriver.PhantomJS(executable_path='/Users/minghao/Research/Softwares/phantomjs-2.1.1-macosx/bin/phantomjs')

#    options = webdriver.ChromeOptions()
#    options.add_argument('headless')
#    options.add_argument('window-size=1920x1080')
#    options.add_argument("disable-gpu")

#    driver = webdriver.Chrome(executable_path='/Users/minghao/Research/Softwares/chromedriver')

    driver.get("http://horus.roe.ac.uk:8080/vdfs/VcrossID_form.jsp")


    if prog=='110':
        database = 'VHSDR6'
    elif prog=='140':
        database = 'VIKINGDR5'


    Select(driver.find_element_by_name("programmeID")).select_by_value(prog)#110 VHS; 140 VIKING
    Select(driver.find_element_by_name("database")).select_by_value(database)
    Select(driver.find_element_by_name("baseTable")).select_by_value('source')#select merged source catalog

    element1 = driver.find_element_by_xpath("//input[@name='selectList']")
    element1.clear()
    parameter = 'SourceID,ra,dec,aJ,aH,aKs,jAperMag3,jAperMag3Err,jAperMagNoAperCorr6,jAperMag6Err,hAperMag3,hAperMag3Err,hAperMagNoAperCorr6,hAperMag6Err,ksAperMag3,ksAperMag3Err,ksAperMagNoAperCorr6,ksAperMag6Err,jErrBits,jppErrBits'
    #parameter = 'SourceID,ra,dec'
    element1.send_keys(parameter)#attributes to select: the parameters

    element2 = driver.find_element_by_xpath("//input[@name='radius']")
    element2.clear()
    element2.send_keys('2.0')

    Select(driver.find_element_by_name("nearest")).select_by_value('1')#1>Nearest object only

    element3 = driver.find_element_by_name("uploadFile")
    element3.send_keys(filename)

    email = 'minghao.astro@gmail.com'
    element4 = driver.find_element_by_name("emailAddress")
#    element4.send_keys(email)

    element5 = driver.find_element_by_xpath("//input[@value='CSV']")#data in ASCII FILE
    element5.click()

    element6 = driver.find_element_by_name("rows")
    element6.clear()
    element6.send_keys('0')

    element7 = driver.find_element_by_xpath("//input[@value='Submit']")
    element7.submit()

    datalink='none'

    for index in range(60):
        driver.implicitly_wait(60)
        #print "after wait 6"
        element_test = driver.find_elements_by_xpath("//a[@href]")
        print("test:",element_test,len(element_test))

        if len(element_test)>0:
            break
        else:
            driver.navigate().refresh()

    for element in driver.find_elements_by_xpath("//a[@href]"):
       #print element.get_attribute('href')[-7:]
        if element.get_attribute('href')[-7:] == '.csv.gz':
            datalink=element.get_attribute('href')
            print(datalink)
            break

    # =======================

    driver.close()
    return datalink

def UHS_CrossID(filename, prog='107'):
    driver = webdriver.PhantomJS(executable_path='/Users/minghao/Research/Softwares/phantomjs-2.1.1-macosx/bin/phantomjs')
#    options = webdriver.ChromeOptions()
#    options.add_argument('headless')
#    options.add_argument('window-size=1920x1080')
#    options.add_argument("disable-gpu")

#    driver = webdriver.Chrome(executable_path='/Users/minghao/Research/Softwares/chromedriver',options=options)
#    driver = webdriver.Chrome(executable_path='/Users/minghao/Research/Softwares/chromedriver')

    driver.get("http://wsa.roe.ac.uk:8080/wsa/crossID_form.jsp")
    driver.refresh()

    if prog=='107':
        database = 'UHSDR1'
    elif prog=='101':
        database = 'UKIDSSDR11PLUS'
    Select(driver.find_element_by_name("database")).select_by_value(database)

    Select(driver.find_element_by_name("programmeID")).select_by_value(prog)#107:UHS; 101:LAS
    Select(driver.find_element_by_name("baseTable")).select_by_value('source')#select merged source catalog

    element1 = driver.find_element_by_xpath("//input[@name='selectList']")
    element1.clear()
    if prog=='101':
        parameter = 'SourceID,ra,dec,ebv,aJ,aH,aK,jAperMag3,jAperMag3Err,jAperMag6,jAperMag6Err,hAperMag3,hAperMag3Err,hAperMag6,hAperMag6Err,kAperMag3,kAperMag3Err,kAperMag6,kAperMag6Err,jErrBits,jppErrBits'
    elif prog=='107':
        parameter = 'SourceID,ra,dec,ebv,aJ,jAperMag3,jAperMag3Err,jAperMag6,jAperMag6Err,jErrBits,jppErrBits'

#    parameter = 'default'
    element1.send_keys(parameter)#attributes to select: the parameters

    element2 = driver.find_element_by_xpath("//input[@name='radius']")
    element2.clear()
    element2.send_keys('2.0')

    Select(driver.find_element_by_name("nearest")).select_by_value('1')#1>Nearest object only

    element3 = driver.find_element_by_name("uploadFile")
    element3.send_keys(filename)

    email = 'minghao.astro@gmail.com'
    element4 = driver.find_element_by_name("emailAddress")
#    element4.send_keys(email)

    element5 = driver.find_element_by_xpath("//input[@value='CSV']")#data in ASCII FILE
    element5.click()

    element6 = driver.find_element_by_name("rows")
    element6.clear()
    element6.send_keys('0')

    element7 = driver.find_element_by_xpath("//input[@value='Submit']")
    element7.click()

    datalink='none'

    for index in range(60):
        print(index)
        driver.implicitly_wait(60)
        #print "after wait 6"
        element_test = driver.find_elements_by_xpath("//a[@href]")
        print("test:",element_test,len(element_test))

        if len(element_test)>0:
            break
#        else:
#            driver.refresh()

    for element in driver.find_elements_by_xpath("//a[@href]"):
       #print element.get_attribute('href')[-7:]
        if element.get_attribute('href')[-7:] == '.csv.gz':
            datalink=element.get_attribute('href')
            print(datalink)
            break

    # =======================

    driver.close()
    return datalink

def download_NIR(tbl, ra_col, dec_col, outdir, survey,\
                 n_onematch=1000):

    if not survey in ['VHS', 'UHS', 'UKIDSS', 'VIKING']:
        print('Survey currently unavailable')
        return 0

    outdir = os.path.abspath(outdir)

    os.system('mkdir -p %s'%outdir)
    os.system('mkdir -p %s/upload'%outdir)
    os.system('mkdir -p %s/download'%outdir)

    Niter = int(len(tbl)/n_onematch)

    for niter in range(Niter+1):
        start = niter * n_onematch
        end = np.min([(niter+1)*n_onematch, len(tbl)])

        output_file = os.path.abspath(outdir+'/download/match%d.csv.gz'%(niter))
        output_file_unzip = output_file[:-3]

        if os.path.exists(output_file_unzip):
            continue

        print('Generating coordinate list for objects %d to %d'%(start, end-1))
        ralist_iter = tbl[ra_col][start:end]
        declist_iter = tbl[dec_col][start:end]

        input_file = os.path.abspath(outdir+'/upload/coord%d.tbl'%(niter))
        np.savetxt(input_file, np.transpose([ralist_iter, declist_iter]), fmt='%.6f')


        if survey=='VHS':
            url = VHS_CrossID(input_file)
        elif survey=='UHS':
            url = UHS_CrossID(input_file)
        elif survey=='UKIDSS':
            url = UHS_CrossID(input_file, prog='101')
        elif survey=='VIKING':
            url = VHS_CrossID(input_file, prog='140')

        os.system('wget -O "%s" "%s"'%(output_file, url))
        os.system('gunzip %s'%output_file)

def match_script(tbl1, fmt1, tbl2, fmt2,\
                ra_col1, dec_col1, ra_col2, dec_col2, outname,\
                suffix_1, suffix_2,\
                icmd1='', icmd2='', ocmd='', radius=2,\
                addcomment=[]):
    # al RA's and DEC's are in deg

    stilts_command = 'java -jar /Users/minghao/Research/Softwares/stilts.jar'

    stilts_command += ' tmatch2'

    stilts_command += ' in1="%s"'%tbl1
    stilts_command += ' in2="%s"'%tbl2

    stilts_command += ' ifmt1="%s"'%fmt1
    stilts_command += ' ifmt2="%s"'%fmt2

    stilts_command += ' matcher="sky"'
    stilts_command += ' params=%d'%radius

    if len(icmd1)>0:
        stilts_command += ' icmd1="%s"'%icmd1
    if len(icmd2)>0:
        stilts_command += ' icmd2="%s"'%icmd2

    stilts_command += ' out=%s'%outname

    stilts_command += ' values1="%s %s"'%(ra_col1, dec_col1)
    stilts_command += ' values2="%s %s"'%(ra_col2, dec_col2)
    stilts_command += ' suffix1="%s"'%(suffix_1)
    stilts_command += ' suffix2="%s"'%(suffix_2)
    stilts_command += ' fixcols=all'

#    stilts_command += ' find=best1'
    for comment in addcomment:
        stilts_command += ' '
        stilts_command += comment

    if len(ocmd)>0:
        stilts_command += ' ocmd="%s"'%ocmd

    return stilts_command

def combine_script(tablelist, formatlist, output, icmd=''):
    stilts_command = 'java -jar /Users/minghao/Research/Softwares/stilts.jar'

    stilts_command += ' tcatn'
    stilts_command += ' nin=%d'%len(tablelist)


    for index in range(len(tablelist)):
        num = index + 1
        stilts_command += ' in%d=%s ifmt%d=%s'%(num, tablelist[index],\
                                               num, formatlist[index])
        if len(icmd)>0:
            stilts_command += ' icmd%d="%s"'%(num, icmd)

    stilts_command += ' out="%s"'%output

    return stilts_command

def main():
    tbl = Table.read('./qso_rsg4p5_v2.0_all_v20180516.fits')
    ra_col = 'radeg'
    dec_col = 'decdeg'

    tbl = tbl[~np.isnan(tbl[ra_col])]

#    download_CatWISE(tbl, ra_col, dec_col, './test_WISE')
#    download_NIR(tbl, ra_col, dec_col, './test_UHS', 'UHS')

    ### filelist
#    direct = os.path.abspath('./test_WISE/download')
#    filenames = os.listdir(direct)

#    fullfilenames = [direct + '/' + fn for fn in filenames]
#    script = combine_script(fullfilenames, ['ipac']*len(fullfilenames), './test_WISE/combined.fits')

#    os.system(script)

    script = match_script(os.path.abspath('./qso_rsg4p5_v2.0_all_v20180516.fits'), 'fits',\
                          os.path.abspath('./test_WISE/combined.fits'), 'fits', \
                          'radeg', 'decdeg', 'ra', 'dec',\
                          os.path.abspath('./matchtest.csv'), '_q','_w')

    os.system(script)

if __name__=='__main__':
    main()
