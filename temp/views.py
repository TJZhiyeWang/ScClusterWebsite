from django.shortcuts import render
from django.shortcuts import HttpResponse
from django.shortcuts import render_to_response
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.csrf import csrf_exempt
from . import Cliq
from . import NMF
from . import BACKSPIN
import os
import time 
import smtplib
from email.mime.text import MIMEText
from email.header import Header
import django.utils.datastructures
from rpy2 import robjects
import zipfile
from dca import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy.api as sc
from dca.api import dca
from scvi.dataset import CsvDataset
from scvi.models import *
from scvi.inference import UnsupervisedTrainer,ClassifierTrainer
import torch
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.externals import joblib
import pandas as pd
import json
from sklearn import metrics
import math
import time
import _thread
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Create your views here.
def help(request):
    return render_to_response("help.html")
@csrf_exempt
def visualupload(request):
    try:
        # 获取上传的文件，如果没有文件，则默认为None  
        dataFile =request.FILES.get("datafile", None)
        labelFile=request.FILES.get("labelfile", None)
        #uniqueNumber前端生成
        uniqueNumber=request.POST['unique']
        #获取当前日期
        today=time.strftime("%Y-%m-%d",time.localtime())
        #为该用户建立一个文件夹,以日期为名方便查找
        if not os.path.exists(os.path.join(BASE_DIR, 'upload/userdata',today)):
            os.mkdir(os.path.join(BASE_DIR, 'upload/userdata',today))
        if not os.path.exists(os.path.join(BASE_DIR, 'upload/userdata',today,uniqueNumber)):
            os.mkdir(os.path.join(BASE_DIR, 'upload/userdata',today,uniqueNumber))
        if not os.path.exists(os.path.join(BASE_DIR, 'upload/userresult',today)):
            os.mkdir(os.path.join(BASE_DIR, 'upload/userresult',today))
        if not os.path.exists(os.path.join(BASE_DIR, 'upload/userresult',today, uniqueNumber)):
            os.mkdir(os.path.join(BASE_DIR, 'upload/userresult',today, uniqueNumber))

        #将取到的数据文件放入该用户的文件夹
        #将取到的文件放入该用户的文件夹
        dataFileName=dataFile.name
        labelFileName=labelFile.name
        dataFilepath=os.path.join(BASE_DIR, 'upload/userdata',today,uniqueNumber,dataFileName)
        lableFilepath=os.path.join(BASE_DIR, 'upload/userdata',today,uniqueNumber,labelFileName)
        destination = open(dataFilepath,'wb+')    # 打开特定的文件进行二进制的写操作  
        for chunk in dataFile.chunks():      # 分块写入文件  
            destination.write(chunk)  
        destination.close() 
        #将取到的标签文件放入该用户的文件夹
        destination = open(lableFilepath,'wb+')    # 打开特定的文件进行二进制的写操作  
        for chunk in labelFile.chunks():      # 分块写入文件  
            destination.write(chunk)  
        destination.close()
        #获取k值 这里mail名字错了，忽略
        kValue=int(request.POST['mail'])
        _thread.start_new_thread(IRGWProcess,(request,dataFilepath,lableFilepath,uniqueNumber,today,kValue,))
    except Exception as e:
        print(e)
        return HttpResponse(json.dumps({"status":"failed"}),content_type="application/json")
    return HttpResponse(json.dumps({"status":"success"}),content_type='application/json')
def IRGWProcess(request,datapath,labelpath,uniqueNumber,today,kValue):
    data=pd.read_csv(datapath,sep=' ')
    clust=pd.read_csv(labelpath,sep=' ')
    tsne = TSNE(n_components=2, learning_rate=100).fit_transform(data.values.T)
    jsondata=[]
    for i in range(min(clust.values[:,1]),max(clust.values[:,1])+1):
        mydict={}
        mydict['name']='group'+str(i)
        mydict['data']=[]
        for j in range(len(clust)):
            if clust.values[j,1]==i:
                xy={}
                xy['xVal']=float(tsne[j,0])
                xy['yVal']=float(tsne[j,1])
                xy['num']=1
                xy['name']=str(clust.values[j,0])
                mydict['data'].append(xy)
        jsondata.append(mydict)
    #存储json结果，list存储为pkl格式,方法为sklearn.externals
    file_root=os.path.join(BASE_DIR,'upload/userresult',today,uniqueNumber,'reduc_dim.pkl')
    joblib.dump(jsondata,file_root)
    del jsondata
    #构造结果
    r = robjects.r
    r.source(os.path.join(BASE_DIR,'upload/inter.R'))
    save_path_root=os.path.join(BASE_DIR,'upload/userresult',today,uniqueNumber)
    rg=r.RG(datapath,kValue,os.path.join(save_path_root,'RGresult.csv'))
    print('generate result file')
    rg=dict(zip(rg.names, map(list,list(rg))))
    r1=metrics.adjusted_rand_score(rg['r1'][1] ,clust.values[:,1])
    r2=metrics.adjusted_rand_score(rg['r2'][1] ,clust.values[:,1])
    r3=metrics.adjusted_rand_score(rg['r3'][1] ,clust.values[:,1])
    r4=metrics.adjusted_rand_score(rg['r4'][1] ,clust.values[:,1])
    r5=metrics.adjusted_rand_score(rg['r5'][1] ,clust.values[:,1])
    ariresult=[r1,r2,r3,r4,r5]
    plt.bar(range(len(ariresult)), ariresult ,color='rgbym',tick_label=['sincera','cidr','sc3','pcareduce','IRGW'])
    plt.ylabel('ARI')
    plt.xlabel('methods')
    #nginx
    #=os.path.join('/var/www/AlgorithmWeb/static/images/user',ip+'ari'+'.png')
    #debug 
    save_path=os.path.join(BASE_DIR,'upload/static/images/user',str(uniqueNumber)+'_bar.png')
    plt.savefig(save_path,dpi = 500)
    sim_M=np.array(rg['sim_M'])
    dim_n=math.sqrt(len(sim_M))
    sim_M=np.reshape(sim_M,(int(dim_n),-1))
    ax = sns.heatmap(sim_M,center=0)
    # save_path=os.path.join('/var/www/AlgorithmWeb/static/images/user',ip+'.png')
    save_path=os.path.join(BASE_DIR,'upload/static/images/user',str(uniqueNumber)+'_heat.png')
    plt.savefig(save_path,dpi = 500)


def IRGWQuery(request):
    try:#用于辨识文件
        uniqueNumber=request.GET['unique']
    except KeyError:
        pass
    #获取当前日期
    today=time.strftime("%Y-%m-%d",time.localtime())
    #only check the last generated file exist
    save_path=os.path.join(BASE_DIR,'upload/static/images/user',str(uniqueNumber)+'_heat.png')
    result={}
    if os.path.exists(save_path):
        ariresult='/static/images/user/'+str(uniqueNumber)+'_bar.png'
        heatdata='/static/images/user/'+str(uniqueNumber)+'_heat.png'
        file_root=os.path.join(BASE_DIR,'upload/userresult',today,uniqueNumber,'reduc_dim.pkl')
        jsondata=joblib.load(file_root)
        result={'status':True,'jsondata':jsondata,'aridata':ariresult,'heatdata':heatdata}
        return HttpResponse(json.dumps(result), content_type="application/json")
    result['status']=False
    return HttpResponse(json.dumps(result),content_type="application/json")

@csrf_exempt
def clustering(request):
    try:
        MailReceiver=request.POST['mail']
    except KeyError:
        pass
    try :
        kValue=int(request.POST['kvalue'])
    except KeyError:
        pass
    try:
        myFile =request.FILES.get("myfile", None)# 获取上传的文件，如果没有文件，则默认为None
        #uniqueNumber前端生成
        uniqueNumber=request.POST['unique']
        #获取当前日期
        today=time.strftime("%Y-%m-%d",time.localtime())
        #为该用户建立一个文件夹,以日期为名方便查找
        if not os.path.exists(os.path.join(BASE_DIR, 'upload/userdata',today)):
            os.mkdir(os.path.join(BASE_DIR, 'upload/userdata',today))
        if not os.path.exists(os.path.join(BASE_DIR, 'upload/userdata',today,uniqueNumber)):
            os.mkdir(os.path.join(BASE_DIR, 'upload/userdata',today,uniqueNumber))
        #将取到的文件放入该用户的文件夹
        filename=myFile.name
        filepath=os.path.join(BASE_DIR, 'upload/userdata',today,uniqueNumber,filename)
        destination = open(filepath,'wb+')    # 打开特定的文件进行二进制的写操作  
        for chunk in myFile.chunks():      # 分块写入文件  
            destination.write(chunk)  
        destination.close()
        _thread.start_new_thread(upload,(request,filepath,filename,uniqueNumber,today,kValue,))
    except Exception as e:
        print(e)
        return HttpResponse(json.dumps({'status':'failed','result':'The backned is busy now, please try later'}),content_type="application/json")
    return HttpResponse(json.dumps({'status':'success','result':'Your data is being processed, please do not close the page'}),content_type="application/json")
def upload(request,filepath,filename,uniqueNumber,today,kValue):
    #首先添加r的算法库，rs为别人得算法，r为IRGW
    rs = robjects.r
    rs.source(os.path.join(BASE_DIR,'upload/AllMethodClust.r'))
    r = robjects.r
    r.source(os.path.join(BASE_DIR,'upload/inter.R'))
    #建立文件夹存储文件
    if not os.path.exists(os.path.join(BASE_DIR, 'upload/userresult',today)):
        os.mkdir(os.path.join(BASE_DIR, 'upload/userresult',today))
    if not os.path.exists(os.path.join(BASE_DIR, 'upload/userresult',today,uniqueNumber)):
        os.mkdir(os.path.join(BASE_DIR, 'upload/userresult',today,uniqueNumber))
    save_path_root=os.path.join(BASE_DIR, 'upload/userresult',today,uniqueNumber)
    #建立一个打包的路径
    zip_root=os.path.join(BASE_DIR,'upload/userresult',today)
    #IRGW
    try:
        ALGOopt=request.POST['cb16']
        if ALGOopt=='on':
            r.RG(filepath,kValue,os.path.join(save_path_root,'IRGW-'+filename))
    except KeyError:
        pass
    #SC3
    try:
        ALGOopt=request.POST['cb6']
        if ALGOopt=='on':
            rs.SC3(filepath,os.path.join(save_path_root,'SC3-'+filename))
            
    except KeyError:
        pass
    #CIDR
    try:
        ALGOopt=request.POST['cb4']
        if ALGOopt=='on':
            rs.CIDR(filepath,os.path.join(save_path_root,'CIDR-'+filename))
            
    except KeyError:
        pass
    #RaceID
    try:
        ALGOopt=request.POST['cb1']
        if ALGOopt=='on':
            rs.RaceID(filepath,os.path.join(save_path_root,'RaceID-'+filename))
            
    except KeyError:
        pass
    #DIMMSC
    try:
        ALGOopt=request.POST['cb13']
        if ALGOopt=='on':
            rs.DIMMSC(filepath,os.path.join(save_path_root,'DIMMSC-'+filename),kValue)
            
    except KeyError:
        pass
    #SINCERA
    try:
        ALGOopt=request.POST['cb5']
        if ALGOopt=='on':
            rs.SINCERA(filepath,os.path.join(save_path_root,'SINCERA-'+filename),kValue)
            
    except KeyError:
        pass
    #Seurat
    try:
        ALGOopt=request.POST['cb9']
        if ALGOopt=='on':
            rs.Seurat(filepath,os.path.join(save_path_root,'Seurat-'+filename))
            
    except KeyError:
        pass
    #cellTree
    try:
        ALGOopt=request.POST['cb3']
        if ALGOopt=='on':
            rs.cellTree(filepath,os.path.join(save_path_root,'cellTree-'+filename),kValue)
            
    except KeyError:
        pass
    #pcaReduce
    try:
        ALGOopt=request.POST['cb2']
        if ALGOopt=='on':
            rs.pcaReduce(filepath,os.path.join(save_path_root,'pcaReduce1-'+filename),kValue,os.path.join(save_path_root,'pcaReduce2-'+filename))
            
    except KeyError:
        pass
    #giniClust
    try:
        ALGOopt=request.POST['cb7']
        if ALGOopt=='on':
            rs.giniClust(filepath,os.path.join(save_path_root,'giniClust-'+filename))
    except KeyError:
        pass
    #SNN
    try:
        ALGOopt=request.POST['cb8']
        if ALGOopt=='on':
            rs.SNN(filepath,os.path.join(save_path_root,'temp-r.txt'))
            Cliq.main(os.path.join(save_path_root,'temp-r.txt'),os.path.join(save_path_root,'temp-py.txt'))
            rs.output(filepath,os.path.join(save_path_root,'temp-py.txt'),os.path.join(save_path_root,'SNN-'+filename))
    except KeyError:
        pass
    #nmf
    try:
        ALGOopt=request.POST['cb10']
        if ALGOopt=='on':
            NMF.FUN_NMF_iterative(filepath,os.path.join(save_path_root,'NMF-'+filename),kValue)
    except KeyError:
        pass
    #BACKSPIN
    try:
        ALGOopt=request.POST['cb11']
        if ALGOopt=='on':
            BACKSPIN.BACKSPIN(filepath,os.path.join(save_path_root,'backspin-'+filename),kValue)
    except KeyError:
        pass
    #biscult
    try:
        ALGOopt=request.POST['cb12']
        if ALGOopt=='on':
            rs.BISCUIT(filepath,os.path.join(save_path_root,'BISCUIT'+filename))
    except KeyError:
        pass
    #dca
    try:
        ALGOopt=request.POST['cb15']
        if ALGOopt=='on':
            data=pd.read_csv(filepath,sep=' ')
            data=data.T
            newfile=os.path.join(BASE_DIR, 'upload/userdata',str(ip),'dca.csv')
            data.to_csv(newfile,sep=',',header=None)
            adata=io.read_dataset(adata=newfile)
            cellname=data.columns.to_list()
            genename=data.index.to_list()
            adata.obs_names = genename
            adata.var_names = cellname
            sc.pp.filter_genes(adata, min_counts=1)
            dca(adata, threads=4)
            y_pred=KMeans(n_clusters=kValue,random_state=9).fit_predict(adata.X)
            result=pd.DataFrame({'name':data.index,'label':y_pred})
            result.to_csv(os.path.join(save_path_root,'DCA'+filename),index=None)
    except KeyError:
        pass
    #scvi
    try:
        ALGOopt=request.POST['cb14']
        if ALGOopt=='on':
            data=pd.read_csv(filepath,sep=' ')
            newfile=os.path.join(BASE_DIR, 'upload/userdata',str(ip),'scvi.csv')
            save_path=os.path.join(BASE_DIR, 'upload/userdata',str(ip))
            data.to_csv(newfile,sep=',',header=None)
            n_epochs_all = None
            show_plot = True
            gene_dataset=CsvDataset(newfile,save_path=save_path,new_n_genes=False)
            n_epochs=40
            lr=1e-3
            use_batches=False
            use_cuda=False
            scvi=SCANVI(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches,n_labels=kValue)
            trainer=UnsupervisedTrainer(scvi,gene_dataset,train_size=1.0)
            trainer.train(n_epochs=n_epochs,lr=lr)
            trainer.posteriors=trainer.create_posterior(indices=np.arange(len(gene_dataset)))
            latent, _, labels = trainer.posteriors.get_latent()
            y_pred=KMeans(n_clusters=kValue,random_state=9).fit_predict(latent)
            result=pd.DataFrame({'name':data.columns,'label':y_pred})
            result.to_csv(os.path.join(save_path_root,'scvi'+filename),index=None)
    except KeyError:
        pass
    #处理完成后发送邮件 todo 带附件发送
    #send_email(MailReceiver)
    #将文件压缩打包
    zip_compress(save_path_root,os.path.join(zip_root,uniqueNumber+'.zip'))
    #返回文件生成成功得消息
    print('the work has been done')
#jquery 轮询判断是否有文件生成，若生成则取回
def isFileExist(request):
    try:#用于辨识文件
        uniqueNumber=request.GET['unique']
    except KeyError:
        pass
    #构造文件路径
    #获取当前日期
    today=time.strftime("%Y-%m-%d",time.localtime())
    zip_root_file=os.path.join(BASE_DIR,'upload/userresult',today,uniqueNumber+'.zip')
    #判断是否有这个文件
    result={}
    if os.path.exists(zip_root_file):
        # targetFile=open(zip_root_file,'rb')
        # response = HttpResponse(targetFile, content_type='application/octet-stream') #设定文件头，这种设定可以让任意文件都能正确下载，而且已知文本文件不是本地打开
        # response['Content-Disposition'] = 'attachment;filename="{0}"'.format(str(uniqueNumber)+'.zip')
        # return response
        result['result']=True
        result['url']=zip_root_file
        return HttpResponse(json.dumps(result),content_type="application/json")
    else:
        result['result']=False
        return HttpResponse(json.dumps(result),content_type="application/json")
        
def downloadFile(request):
    filename=request.GET['url']
    targetFile=open(filename,'rb')
    response=HttpResponse(targetFile)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format('result.zip')  
    return response

def visual(request):
    return render(request, 'methods.html')
def index(request):
    return render(request, 'index.html')
    
def datasets(request):
    #DataFileName=file_name(os.path.join(BASE_DIR, 'upload/datasets/data'))
    #file_name('C:/Users/49387/AlgorithmWeb/upload/datasets/data')
    #ClustFileName=file_name(os.path.join(BASE_DIR, 'upload/datasets/clust'))
    #file_name('C:/Users/49387/AlgorithmWeb/upload/datasets/clust')
    #FeedBack={'DataFileName':DataFileName,'ClustFileName':ClustFileName}
    return render_to_response("datasets.html")
    
def download(request,name):
    file=open(os.path.join(BASE_DIR,'upload/datasets',name+'.zip'),'rb')
    #open('C:/Users/49387/AlgorithmWeb/upload/datasets/'+name,'rb')
    response=HttpResponse(file)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(name+'.zip')  
    return response 
@csrf_protect
def cluster(request):
   
    return render_to_response("clustering.html")
def process(destination):
    global rateOfprocess
    while rateOfprocess < 10 :
        time.sleep(1)
        rateOfprocess+=1
    rateOfprocess=0

def show_progress(request):
    global rateOfprocess
    return JsonResponse(rateOfprocess, safe=False)
    
    
    #给指定的邮箱发送一封邮件
def send_email(receiver):
    username='ant496@163.com'
    password='cell123456'
    message = MIMEText('CONGRATULATIONS! YOUR DATA HAS BEEN PROCESSED SUCCESSFULLY. PLEASE DOWNLOAD ON OUR WEBSITE IN TIME!', 'plain', 'utf-8')
    message['From'] = 'TJCS<ant496@163.com>'   # 发送者
    message['To'] = 'USER<'+ receiver +'>'    # 接收者
    subject = 'Your data has been processed successfully!!!'
    message['Subject'] = Header(subject, 'utf-8')
    
    
    smtp = smtplib.SMTP()
    smtp.connect('smtp.163.com',25)
    
    smtp.login(username, password) 
    smtp.sendmail(username, receiver, message.as_string()) 
    smtp.quit()
    #获取某个路径下的所有的文件名称返回列表
def file_name(file_dir):
    L=[]   
    for root, dirs, files in os.walk(file_dir):  
        for file in files:  
            L.append(file)  
    return L  
#压缩文件
def zip_compress(startdir,filenews):
    z = zipfile.ZipFile(filenews,'w',zipfile.ZIP_DEFLATED) #参数一：文件夹名
    for dirpath, dirnames, filenames in os.walk(startdir):
        for filename in filenames:
            z.write(os.path.join(dirpath,filename),filename)
    z.close()
