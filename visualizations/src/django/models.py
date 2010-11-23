from django.db import models

class CancerType(models.Model):
    name = models.CharField(max_length=30)

class Platform(models.Model):
    name = models.CharField(max_length=50)
    probeset_path = models.CharField(max_length=50)

class ExpressionData(models.Model):
    path = models.CharField(max_length=50)
    platform = models.OneToOneField(Platform)

class SurvivalData(models.Model):
    path = models.CharField(max_length=50)
    platform = models.OneToOneField(Platform)

class DataSource(models.Model):
    name = models.CharField(max_length=30)
    description = models.CharField(max_length=500)
    menu_text = models.CharField(max_length=15)
    
