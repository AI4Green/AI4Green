@description('The name of the container app environment')
param containerAppEnvName string

@description('The location of the container app environment')
param location string = resourceGroup().location

param tags object = {}

// Azure Container App Environment (required for container apps)
resource containerEnv 'Microsoft.App/managedEnvironments@2023-05-01' = {
  name: containerAppEnvName
  location: location
  properties: {}
  tags: union({
    Source: 'Bicep'
  }, tags)
}

output id string = containerEnv.id
