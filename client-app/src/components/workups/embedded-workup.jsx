import React, { useState } from "react";

import { WorkupForm } from "components/workups/form";
import { WorkupCreateInline } from "components/workups/create";
import { WorkupCanvas } from "pages/workups/canvas";

const WORKUP_VIEWS = {
  SELECT_TEMPLATE: "select-template",
  CREATE_TEMPLATE: "create-template",
  EDIT_TEMPLATE: "edit-template",
};

export const EmbeddedWorkup = ({
  reactionId,
  initialInstanceId = null,
  workbookName,
  workgroupName,
}) => {
  const [instanceId, setInstanceId] = useState(initialInstanceId);
  const [view, setView] = useState("select");

  if (instanceId) {
    return <WorkupForm instanceId={instanceId} />;
  }

  if (view === "template") {
    return (
      <WorkupCanvas
        onCancel={() => setView("select")}
        onCreated={(template) => {
          // Decide whether you want to immediately instantiate
          // it or return to template selection.
          setView("select");
        }}
      />
    );
  }

  return (
    <WorkupCreateInline
      reactionId={reactionId}
      workgroupName={workgroupName}
      workbookName={workbookName}
      onCreateTemplate={() => setView("template")}
      onCreated={(data) => {
        setInstanceId(data.id);
      }}
    />
  );
};
