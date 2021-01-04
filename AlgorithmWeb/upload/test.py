import smtplib
from email.mime.text import MIMEText
from email.header import Header
import numpy as np
def send_email(receiver):
    username='ant496@163.com'
    password='cell123456'
    message = MIMEText('hello,my name is wangzhiye, i am from tongji university and this is a test-email.', 'plain', 'utf-8')
    message['From'] = 'tongji<ant496@163.com>'   # 发送者
    message['To'] = 'user<'+ receiver +'>'    # 接收者
    subject = 'Your data has been processed successfully!!!'
    message['Subject'] = Header(subject, 'utf-8')
    
    
    smtp = smtplib.SMTP()
    smtp.connect('smtp.163.com',25)
    
    smtp.login(username, password) 
    smtp.sendmail(username, receiver, message.as_string()) 
    smtp.quit()
    '''try:
        smtp = smtplib.SMTP()
        smtp.connect('smtp.163.com',25)
        #smtp = smtplib.SMTP_SSL('smtp.163.com', port=465)
        smtp.login(username, password) 
        smtp.sendmail(username, receiver, message.as_string()) 
        smtp.quit()
    except smtplib.SMTPException:
        print ("Error: 无法发送邮件")'''
#send_email('493878030@qq.com')
import os  
def file_name(file_dir):
    L=[]   
    for root, dirs, files in os.walk(file_dir):  
        for file in files:  
            L.append(file)  
    return L  
s=np.array([[  1,   2,   3, 100],
       [  4,   5,   6, 200],
       [  7,   8,   9, 300]])
print(s[1:,1:])