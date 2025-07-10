// adapted from https://github.com/uon-drs/bicep/blob/main/modules/components/vnet.bicep

param vnetName string

param location string = resourceGroup().location

param tags object = {}


resource vnet 'Microsoft.Network/virtualNetworks@2020-08-01' = {
  name: vnetName
  location: location
  properties: {
    addressSpace: {
      addressPrefixes: [
        '10.0.0.0/23'
      ]
    }
    // Microsoft.Network/virtualNetworks/subnets
    subnets: [
      {
        name: 'Default' // for VM, most clients
        properties: {
          addressPrefix: '10.0.0.0/23'
          privateEndpointNetworkPolicies: 'Enabled'
          privateLinkServiceNetworkPolicies: 'Enabled'
        }
      }
      // {
      //   name: 'Integration' // for App Service VNET Integrations
      //   properties: {
      //     addressPrefix: '10.0.0.0/23'
      //     privateEndpointNetworkPolicies: 'Enabled'
      //     privateLinkServiceNetworkPolicies: 'Enabled'
      //     delegations: [
      //       {
      //         name: 'delegation'
      //         properties:{
      //           serviceName: 'Microsoft.Web/serverFarms'
      //         }
      //       }
      //     ]
      //   }
      // }
    ]
  }
  tags: union({
    Source: 'Bicep'
  }, tags)
}

output name string = vnet.name
output defaultSubnetId string = vnet.properties.subnets[0].id
// output integrationSubnetId string = vnet.properties.subnets[1].id
