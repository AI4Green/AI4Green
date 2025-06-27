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
          env: [
            {
              name: 'KAFKA_BROKER_ID'
              value: '1'
            }
            {
              name: 'KAFKA_ZOOKEEPER_CONNECT'
              value: 'TODO' // get url for the zooper image
            }
            {
              name: 'KAFKA_LISTENER_SECURITY_PROTOCOL_MAP'
              value: 'PLAINTEXT:PLAINTEXT,PLAINTEXT_HOST:PLAINTEXT'
            }
            {
              name: 'KAFKA_ADVERTISED_LISTENERS'
              // TODO: substitue the kafka url with the url to this container in azure
              value: 'PLAINTEXT://kafka:29092,PLAINTEXT_HOST://localhost:9092'
            }
            {
              name: 'KAFKA_LISTENERS'
              value: 'PLAINTEXT://0.0.0.0:29092,PLAINTEXT_HOST://0.0.0.0:9092'
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
          name: '${appBaseNae}-kafka-connect'
          image: kafkaConnectImage
        }
      ]
    }
  }
}
