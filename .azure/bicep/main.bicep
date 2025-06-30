import { referenceSecret } from 'br/DrsUtils:functions:v1'
import { ConnectionStringDictionary } from 'br/DrsUtils:types:v1'

type ServiceNames = 'aigreen'
param serviceName ServiceNames

type Environments = 'dev' | 'qa' | 'uat' | 'prod'
param env Environments

param keyVaultName string = env == 'uat' ? '${serviceName}-test-kv' : '${serviceName}-${env}-kv'

param appName string = '${serviceName}-${env}'
param appSettings object = {}
param hostnames array = []

param location string = resourceGroup().location

param sharedEnv string = 'shared'
var sharedPrefix = serviceName

param appServicePlanSku string = 'P2'


// MinIO params
param minioImage string

// Kafka stack params
param kafkaConnectImage string

// log analytics workspace
module la 'br/DrsComponents:log-analytics-workspace:v1' = {
  name: 'la-ws-${uniqueString(sharedPrefix)}'
  params: {
    location: location
    logAnalyticsWorkspaceName: '${sharedPrefix}-la-ws'
    tags: {
      ServiceScope: serviceName
      Environment: sharedEnv
    }
  }
}

// Provision virtual network
module vnet 'br/DrsComponents:vnet:v1' = {
  params: {
    vnetName: '${serviceName}-${env}-vnet'
    location: location
    tags: {
      ServiceScope: serviceName
      Environment: sharedEnv
    }
  }
}

// keyvault
resource kv 'Microsoft.KeyVault/vaults@2019-09-01' existing = {
  name: keyVaultName
}

// App Service Plan
module asp 'br/DrsComponents:app-service-plan:v1' = {
  name: 'asp'
  params: {
    location: location
    aspName: 'ai4green-asp'
    sku: appServicePlanSku
    tags: {
      ServiceScope: serviceName
      Environment: sharedEnv
    }
  }
}


// App service
module ai4greenSite 'br/DrsComponents:app-service:v1' = {
  params: {
    appName: appName
    aspName: asp.outputs.name
    logAnalyticsWorkspaceName: la.outputs.name
    appHostnames: hostnames
    location: location
  }
}

// Site config
module ai4greenSiteConfig 'br/DrsConfig:webapp:v1' = {
  name: 'siteConfig-${uniqueString(ai4greenSite.name)}'
  params: {
    appName: ai4greenSite.outputs.name
    appFramework: 'PYTHON|3.12'
    startCommand: 'startup.sh' // nextjs standalone
    appSettings: union(appInsightsSettings, ai4greenEnvVars, appSettings)
  }
}

var appInsightsSettings = {
  // App Insights
  ApplicationInsightsAgent_EXTENSION_VERSION: '~2'
  XDT_MicrosoftApplicationInsights_Mode: 'recommended'
  DiagnosticServices_EXTENSION_VERSION: '~3'
  APPINSIGHTS_PROFILERFEATURE_VERSION: '1.0.0'
  APPINSIGHTS_SNAPSHOTFEATURE_VERSION: '1.0.0'
  InstrumentationEngine_EXTENSION_VERSION: '~1'
  SnapshotDebugger_EXTENSION_VERSION: '~1'
  XDT_MicrosoftApplicationInsights_BaseExtensions: '~1'
}

// AI4Green web app environment settings 
var ai4greenEnvVars = {
  APPLICATIONINSIGHTS_CONNECTION_STRING: ai4greenSite.outputs.appInsights.connectionString
  OIDC_CLIENT_ID: referenceSecret(keyVaultName, 'OIDC_CLIENT_ID')
  OIDC_CLIENT_SECRET: referenceSecret(keyVaultName, 'OIDC_CLIENT_SECRET')
  OIDC_CLIENT_AUTH_URI: referenceSecret(keyVaultName, 'OIDC_CLIENT_AUTH_URI')
  OIDC_CLIENT_TOKEN_URI: referenceSecret(keyVaultName, 'OIDC_CLIENT_TOKEN_URI')
  OIDC_USERINFO_AUTH_URI: referenceSecret(keyVaultName, 'OIDC_CLIENT_USERINFO_URI')
  OIDC_CLIENT_ISSUER: referenceSecret(keyVaultName, 'OIDC_CLIENT_ISSUER')
  OIDC_CLIENT_REDIRECT_URIS: referenceSecret(keyVaultName, 'OIDC_CLIENT_REDIRECT_URIS')
  MESSAGE_QUEUE_HOSTNAME: kafkaStack.outputs.kafkaUrl
  USE_KAFKA: 1
  MINIO_ACCESS_KEY: minio.outputs.accessKey
  MINIO_SECRET_KEY: minio.outputs.secretKey
  MINIO_SECURE: 1
}

// Grant frontend Key Vault access
module ai4greendKvAccess 'br/DrsConfig:keyvault-access:v2' = {
  name: 'kvAccess-${uniqueString(ai4greenSite.name)}'
  params: {
    keyVaultName: kv.name
    principalId: ai4greenSite.outputs.identity.principalId
  }
}

// Provision container app managed environment
module containerAppEnv 'container-env.bicep' = {
  params: {
    containerAppEnvName: '${serviceName}-${env}-env'
    location: location
    infrastructureSubnetId: vnet.outputs.integrationSubnetId
  }
}


// Provision Minio
module minio 'minio.bicep' = {
  params: {
    location: location
    storageAccountName: '${serviceName}-${env}-storage-account'
    fileShareName: '${serviceName}-${env}-file-share'
    containerImage: minioImage
    containerAppName: '${serviceName}-${env}-minio'
    containerAppEnvId: containerAppEnv.outputs.id
    minioAccessKey: referenceSecret(keyVaultName, 'MINIO_CLIENT_ACCESS_KEY')
    minioBucketName: referenceSecret(keyVaultName, 'MINIO_BUCKET_NAME')
    minioRootPassword: referenceSecret(keyVaultName, 'MINIO_ROOT_PASSWORD')
    minioRootUserName: referenceSecret(keyVaultName, 'MINIO_ROOT_USER')
    minioSecretKey: referenceSecret(keyVaultName, 'MINIO_CLIENT_SECRET_KEY')
  }
}


// Provision kafka stack
module kafkaStack 'kafka.bicep' = {
  params: {
    appBaseNae: '${serviceName}-${env}'
    location: location
    kafkaConnectImage: kafkaConnectImage
    containerAppEnvId: containerAppEnv.outputs.id
    minioAccessKey: minio.outputs.accessKey
    minioSecretKey: minio.outputs.secretKey
    minioUrl: minio.outputs.appUrl
  }
}
