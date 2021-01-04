from django.contrib import admin
from django.urls import path
from upload import views
from django.urls import re_path
app_name = 'upload'

urlpatterns = [
    re_path('admin/', admin.site.urls),
    re_path('index/',views.index),
    re_path('datasets/',views.datasets),
	re_path(r'^download/(?P<name>\S+)/$',views.download,name='download'),
    re_path('clustering/',views.cluster),
    re_path(r'^upload/$',views.upload),
    re_path('visulization/',views.visual),
    re_path('visualupload/',views.visualupload),
    re_path('help/',views.help),
]
