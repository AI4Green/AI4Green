@description('Location for all resources.')
param location string = resourceGroup().location

@description('Base name for the apps')
param appBaseNae string

@description('The image to use for zookeeper')
param zookeeperImage string = 'confluentinc/cp-zookeeper:7.6.0'

@description('The image to use for kafka')
param kafkaImage string = 'confluentinc/cp-kafka:7.6.0'

@description('The image to use for kafka-connect')
param kafkaConnectImage string = 'confluentinc/cp-kafka-connect:7.6.0'

// Zookeeper settings
var zookeeperClientPort string = '2181'
var zookeeperTickTime string = '2000'

// Kafka settings


// Kafka Connect settings


// Kafka stack container app
resource containerApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: '${appBaseNae}-kafka-stack'
  location: location
  properties: {
    template: {
      containers: [
        // Zookeeper
        {
          name: '${appBaseNae}-zookeeper'
          image: zookeeperImage
          env: [
            {
              name: 'ZOOKEEPER_CLIENT_PORT'
              value: zookeeperClientPort
            }
            {
              name: 'ZOOKEEPER_TICK_TIME'
              value: zookeeperTickTime
            }
          ]
        }
        // Kafka
        {
          name: '${appBaseNae}-kafka'
          image: kafkaImage
        }
        // Kafka-connect
        {
          name: '${appBaseNae}-kafka-connect'
          image: kafkaConnectImage
        }
      ]
    }
  }
}
