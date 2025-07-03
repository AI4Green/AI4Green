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

// Zookeeper settings
var zookeeperClientPort int = 2181
var zookeeperTickTime int = 2000

// Kafka settings
var kafkaMainPort int = 9092
var kafkaAltPort int = 29092


// Kafka Connect settings

// Kafka stack container app
resource containerApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: '${appBaseName}-kafka-stack'
  location: location
  properties: {
    managedEnvironmentId: containerAppEnvId
    configuration: {
      ingress: {
        exposedPort: kafkaMainPort
        targetPort: kafkaMainPort
        external: true
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
              value: 'PLAINTEXT://localhost:${kafkaAltPort},PLAINTEXT_HOST://localhost:${kafkaMainPort}'
            }
            {
              name: 'KAFKA_LISTENERS'
              value: 'PLAINTEXT://0.0.0.0:${kafkaAltPort},PLAINTEXT_HOST://0.0.0.0:${kafkaMainPort}'
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
        // Kafka-connect
        {
          name: '${appBaseName}-kafka-connect'
          image: kafkaConnectImage
          env: [
            {
              name: 'CONNECT_BOOTSTRAP_SERVERS'
              value: 'localhost:${kafkaAltPort}'
            }
            {
              name: 'BOOTSTRAP_SERVERS'
              value: 'localhost:${kafkaAltPort}'
            }
          ]
        }
      ]
    }
  }
}

output kafkaUrl string = 'https://${containerApp.properties.configuration.ingress.fqdn}:${kafkaMainPort}'
