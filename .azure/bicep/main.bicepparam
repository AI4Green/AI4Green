using './main.bicep'

param serviceName = 'ai4green'
param env = 'uat'
param minioImage = 'ai4greeneln/ai4green-minio'
param kafkaConnectImage = 'ai4greeneln/ai4green-kafka-connect'
param hostnames = ['ai4green-uat.azurewebsites.net']
