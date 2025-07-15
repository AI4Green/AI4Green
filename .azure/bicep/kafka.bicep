@description('Location for all resources.')
param location string = resourceGroup().location

@description('Base name for the apps')
param appBaseName string

@description('The image to use for zookeeper')
param zookeeperImage string = 'confluentinc/cp-zookeeper:7.6.0'

@description('The image to use for kafka')
param kafkaImage string = 'confluentinc/cp-kafka:7.6.0'

@description('The image to use for kafka-connect')
param kafkaConnectImage string = 'confluentinc/cp-kafka-connect:7.6.0'

@description('Resource ID of the managed environment for Azure Container Apps.')
param containerAppEnvId string

@description('Name of the managed environment for Azure Container Apps.')
param containerAppEnvName string

param tags object = {}

// Zookeeper settings
var zookeeperClientPort int = 2181
var zookeeperTickTime int = 2000

// Kafka settings
var kafkaMainPort int = 9092
var kafkaAppName string = '${appBaseName}-kafka'

// Kafka Connect settings
var kafkaConnectPort int = 8083
var kafkaConnectAppName string = '${appBaseName}-kafka-connect'

// Kafka container app
resource kafkaApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: kafkaAppName
  location: location
  properties: {
    managedEnvironmentId: containerAppEnvId
    configuration: {
      ingress: {
        targetPort: kafkaMainPort
      }
    }
    template: {
      containers: [
        // Zookeeper
        {
          name: '${appBaseName}-zookeeper'
          image: zookeeperImage
          env: [
            {
              name: 'ZOOKEEPER_CLIENT_PORT'
              value: string(zookeeperClientPort)
            }
            {
              name: 'ZOOKEEPER_TICK_TIME'
              value: string(zookeeperTickTime)
            }
          ]
        }
        // Kafka
        {
          name: '${appBaseName}-kafka'
          image: kafkaImage
          env: [
            {
              name: 'KAFKA_BROKER_ID'
              value: '1'
            }
            {
              name: 'KAFKA_ZOOKEEPER_CONNECT'
              value: 'localhost:${zookeeperClientPort}'
            }
            {
              name: 'KAFKA_LISTENER_SECURITY_PROTOCOL_MAP'
              value: 'PLAINTEXT:PLAINTEXT,PLAINTEXT_HOST:PLAINTEXT'
            }
            {
              name: 'KAFKA_ADVERTISED_LISTENERS'
              value: 'PLAINTEXT_HOST://${kafkaAppName}.${containerAppEnvName}.internal:${kafkaMainPort}'
            }
            {
              name: 'KAFKA_LISTENERS'
              value: 'PLAINTEXT_HOST://0.0.0.0:${kafkaMainPort}'
            }
            {
              name: 'KAFKA_INTER_BROKER_LISTENER_NAME'
              value: 'PLAINTEXT'
            }
            {
              name: 'KAFKA_OFFSETS_TOPIC_REPLICATION_FACTOR'
              value: '1'
            }
          ]
        }
      ]
    }
  }
  tags: union({
    Source: 'Bicep'
  }, tags)
}

// Kafka Connect container app
resource kafkaConnectApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: kafkaConnectAppName
  location: location
  properties: {
    managedEnvironmentId: containerAppEnvId
    configuration: {
      ingress: {
        targetPort: kafkaConnectPort
      }
    }
    template: {
      containers: [
        // Zookeeper
        {
          name: '${appBaseName}-zookeeper'
          image: zookeeperImage
          env: [
            {
              name: 'ZOOKEEPER_CLIENT_PORT'
              value: string(zookeeperClientPort)
            }
            {
              name: 'ZOOKEEPER_TICK_TIME'
              value: string(zookeeperTickTime)
            }
          ]
        }
        // Kafka-connect
        {
          name: '${appBaseName}-kafka-connect'
          image: kafkaConnectImage
          env: [
            {
              name: 'CONNECT_BOOTSTRAP_SERVERS'
              value: '${kafkaApp.name}.${containerAppEnvName}.internal:${kafkaMainPort}'
            }
            {
              name: 'BOOTSTRAP_SERVERS'
              value: '${kafkaApp.name}.${containerAppEnvName}.internal:${kafkaMainPort}'
            }
          ]
        }
      ]
    }
  }
  tags: union({
    Source: 'Bicep'
  }, tags)
}

output kafkaUrl string = 'https://${kafkaApp.properties.configuration.ingress.fqdn}:${kafkaMainPort}'
output kafkaConnectUrl string = 'https://${kafkaConnectApp.properties.configuration.ingress.fqdn}:${kafkaMainPort}'
