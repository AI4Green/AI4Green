@description('The name of the container app environment')
param containerAppEnvName string

@description('The location of the container app environment')
param location string = resourceGroup().location

@description('The ID of the virtual network')
param infrastructureSubnetId string

// Azure Container App Environment (required for container apps)
resource containerEnv 'Microsoft.App/managedEnvironments@2023-05-01' = {
  name: containerAppEnvName
  location: location
  properties: {
    vnetConfiguration: {
      infrastructureSubnetId: infrastructureSubnetId
    }
  }
}

output id string = containerEnv.id
