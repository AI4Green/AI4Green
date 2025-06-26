@description('Location for all resources.')
param location string = resourceGroup().location

@description('Name of the Container App.')
param containerAppName string = 'minio-app'

@description('MinIO container image to deploy.')
param containerImage string = 'quay.io/minio/minio:RELEASE.2025-04-22T22-12-26Z'

@description('Globally unique name of the Azure Storage Account.')
param storageAccountName string

@description('Name of the Azure File Share to attach to the container.')
param fileShareName string = 'minioshare'

@description('Name of the managed environment for Azure Container Apps.')
param containerAppEnvName string = 'minio-env'

@description('The username of the MinIO root user')
param minioRootUserName string

@description('The password of the MinIO root user')
@secure()
param minioRootPassword string

@description('The name of the bucket to create')
param minioBucketName string

@description('The access key to the bucket')
param minioAccessKey string

@description('The secret key to the bucket')
#disable-next-line secure-secrets-in-params
param minioSecretKey string

var storageMountPath = '/data' // Mount point inside the container
var minioPort = 9000

// Storage account for persistent data
module storageAccount 'br/DrsComponents:storage-account:v1' = {
  params: {
    baseAccountName: storageAccountName
  }
}

// Azure File Share inside the storage account
resource fileShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2021-09-01' = {
  name: '${storageAccount.name}/default/${fileShareName}'
  properties: {
    shareQuota: 100 // Max quota in GB
  }
}

// User-assigned managed identity to allow the container app to access resources securely
resource userAssignedIdentity 'Microsoft.ManagedIdentity/userAssignedIdentities@2023-01-31' = {
  name: '${containerAppName}-identity'
  location: location
}

// Azure Container App Environment (required for container apps)
resource containerEnv 'Microsoft.App/managedEnvironments@2023-05-01' = {
  name: containerAppEnvName
  location: location
  properties: {
    daprAIInstrumentationKey: '' // Optional, set if using Dapr with App Insights
  }
}

// The main MinIO container app with Azure File volume mount
resource containerApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: containerAppName
  location: location
  identity: {
    type: 'UserAssigned'
    userAssignedIdentities: {
      '${userAssignedIdentity.id}': {} // Attach identity
    }
  }
  properties: {
    managedEnvironmentId: containerEnv.id
    configuration: {
      activeRevisionsMode: 'Single'
      secrets: [
        {
          name: 'storage-account-name'
          value: storageAccount.name
        }
      ]
    }
    template: {
      containers: [
        {
          name: 'minio'
          image: containerImage
          env: [
            {
              name: 'MINIO_ROOT_USER'
              value: minioRootUserName
            }
            {
              name: 'MINIO_ROOT_PASSWORD'
              value: minioRootPassword
            }
            {
              name: 'MINIO_BUCKET_NAME'
              value: minioBucketName
            }
            {
              name: 'MINIO_CLIENT_ACCESS_KEY'
              value: minioAccessKey
            }
            {
              name: 'MINIO_CLIENT_SECRET_KEY'
              value: minioSecretKey
            }
            {
              name: 'MINIO_STORAGE_MOUNT'
              value: storageMountPath
            }
          ]
          volumeMounts: [
            {
              volumeName: 'miniodata'
              mountPath: storageMountPath
            }
          ]
        }
      ]
      volumes: [
        {
          name: 'miniodata'
          storageType: 'EmptyDir'
        }
      ]
    }
  }
}

// Export the public FQDN of the MinIO Container App
output appUrl string = 'https://${containerApp.properties.configuration.ingress.fqdn}:${minioPort}'
// Export the access key and secret key for other apps to use
output accessKey string = minioAccessKey
output secretKey string = minioSecretKey
