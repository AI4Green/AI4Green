import "github-markdown-css";

import {
  Container,
  Heading,
  ListItem,
  OrderedList,
  UnorderedList,
  VStack,
} from "@chakra-ui/react";
import { AutoLink } from "components/core/auto-link";
import ky from "ky";
import ReactMarkdown from "react-markdown";
import rehypeRaw from "rehype-raw";
import gfm from "remark-gfm";
import useSWR from "swr";

export const ContentPage = ({ contentKey }) => {
  const { data: content } = useSWR(
    `/locales/dev/${contentKey}.md`,
    async (url) => ky.get(url).text(),
    {
      suspense: true,
    },
  );

  return (
    <Container maxW="75%">
      <VStack
        align="start"
        p={2}
        spacing={4}
        mt={8}
        mb={4}
        fontSize="sm"
        fontWeight="thin"
        letterSpacing="wide"
      >
        <ReactMarkdown
          remarkPlugins={[gfm]}
          rehypePlugins={[rehypeRaw]}
          components={{
            // replace some plain HTML with Chakra components
            // which nets us desirable styling mostly
            a: ({ href, ...props }) => <AutoLink url={href} {...props} />,
            h1: (props) => <Heading size="xl" {...props} />,
            h2: (props) => <Heading size="lg" {...props} />,
            h3: (props) => <Heading size="md" {...props} />,
            h4: (props) => <Heading size="md" {...props} />,
            h5: (props) => <Heading size="sm" {...props} />,
            h6: (props) => <Heading size="xs" {...props} />,
            ul: (props) => <UnorderedList pl={8} {...props} />,
            ol: (props) => <OrderedList pl={8} {...props} />,
            li: (props) => <ListItem {...props} />,
          }}
        >
          {content}
        </ReactMarkdown>
      </VStack>
    </Container>
  );
};
