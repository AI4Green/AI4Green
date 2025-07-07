FROM confluentinc/cp-kafka-connect:7.6.0

# Default envrionment variables
ENV CONNECT_PLUGIN_PATH=/usr/share/java,/etc/kafka-connect/jars
ENV CONNECT_REST_PORT=8083
ENV CONNECT_GROUP_ID=kafka-connect-group
ENV CONNECT_CONFIG_STORAGE_TOPIC=connect-configs
ENV CONNECT_OFFSET_STORAGE_TOPIC=connect-offsets
ENV CONNECT_STATUS_STORAGE_TOPIC=connect-status
ENV CONNECT_KEY_CONVERTER=org.apache.kafka.connect.storage.StringConverter
ENV CONNECT_KEY_CONVERTER_SCHEMAS_ENABLE=false
ENV CONNECT_VALUE_CONVERTER=org.apache.kafka.connect.json.JsonConverter
ENV CONNECT_VALUE_CONVERTER_SCHEMAS_ENABLE=false
ENV CONNECT_INTERNAL_KEY_CONVERTER=org.apache.kafka.connect.json.JsonConverter
ENV CONNECT_INTERNAL_VALUE_CONVERTER=org.apache.kafka.connect.json.JsonConverter
ENV CONNECT_LOG4J_LOGGERS=org.reflections=ERROR
ENV CONNECT_REST_ADVERTISED_HOST_NAME=localhost

USER root

# Install necessary utilities
RUN yum install -y curl jq unzip

# Create plugin directory if it doesn't exist
RUN mkdir -p /etc/kafka-connect/jars

# Copy the connector zip file and unzip it into the plugin path
COPY confluentinc-kafka-connect-s3-10.6.7.zip /tmp/s3-sink-connector.zip
RUN unzip /tmp/s3-sink-connector.zip -d /etc/kafka-connect/jars && rm /tmp/s3-sink-connector.zip

# Copy the custom entrypoint script
COPY kafka-connect-setup.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Kafka Connect will start as normal with the S3 sink plugin loaded