"""AlgorithmWeb URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from upload import views
import upload.urls
from django.urls import re_path
from django.conf.urls import url,include
from AlgorithmWeb import settings
urlpatterns = [
    url('upload/',include(upload.urls,namespace='upload')),
    ]
'''path('admin/', admin.site.urls),
    path('index/',views.index),
    path('datasets/',views.datasets,name='upload'),
	path(r'^download/$)',views.download),
	path('wait/',views.show_progress),'''
