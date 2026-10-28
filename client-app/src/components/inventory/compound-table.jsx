import { Fragment, useState } from "react";
import ContainersTable, { Badge } from "./containers-table.jsx";

const COLUMN_COUNT = 8;

/**
 * rows: [{ wc, summary, items }]
 *   wc       the WorkgroupCompound
 *   summary  result of summarise(wc)
 *   items    containers to show when the row is expanded (already filtered)
 */
export default function CompoundTable({ rows, onItemAction, onAddContainer }) {
  const [expanded, setExpanded] = useState(() => new Set());

  const toggle = (id) =>
    setExpanded((prev) => {
      const next = new Set(prev);
      next.has(id) ? next.delete(id) : next.add(id);
      return next;
    });

  return (
    <table className="inv-table">
      <thead>
        <tr>
          <th className="inv-col-toggle">
            <span className="inv-sr">Expand</span>
          </th>
          <th>Compound</th>
          <th>Hazards</th>
          <th>Total stock</th>
          <th>Containers</th>
          <th>Locations</th>
          <th>Status</th>
          <th>
            <span className="inv-sr">Actions</span>
          </th>
        </tr>
      </thead>
      <tbody>
        {rows.map(({ wc, summary, items }) => {
          const isOpen = expanded.has(wc.id);
          const name = wc.compound.name;

          return (
            <Fragment key={wc.id}>
              <tr className={isOpen ? "inv-row-open" : undefined}>
                <td>
                  <button
                    type="button"
                    className="inv-expand"
                    aria-expanded={isOpen}
                    aria-label={`${isOpen ? "Collapse" : "Expand"} ${name}`}
                    onClick={() => toggle(wc.id)}
                  >
                    <svg
                      width="16"
                      height="16"
                      viewBox="0 0 16 16"
                      aria-hidden="true"
                    >
                      <polyline
                        points={isOpen ? "3,6 8,11 13,6" : "6,3 11,8 6,13"}
                        fill="none"
                        stroke="currentColor"
                        strokeWidth="2"
                      />
                    </svg>
                  </button>
                </td>
                <td>
                  <b>{name}</b>
                  <div className="inv-cas">{wc.compound.cas}</div>
                </td>
                <td>
                  {wc.compound.hazards.map((h) => (
                    <span key={h} className="inv-hazard">
                      {h}
                    </span>
                  ))}
                </td>
                <td>{summary.totalLabel}</td>
                <td>{summary.containerCount}</td>
                <td>{summary.locations.join(", ") || "-"}</td>
                <td>
                  <Badge tone={summary.status.tone}>
                    {summary.status.label}
                  </Badge>
                </td>
                <td>
                  <button
                    type="button"
                    className="inv-link"
                    aria-label={`Add container of ${name}`}
                    onClick={() => onAddContainer(wc)}
                  >
                    Add container
                  </button>
                </td>
              </tr>
              {isOpen && (
                <tr>
                  <td colSpan={COLUMN_COUNT} className="inv-expansion">
                    <div className="inv-sub">
                      <ContainersTable
                        items={items}
                        onItemAction={onItemAction}
                      />
                    </div>
                  </td>
                </tr>
              )}
            </Fragment>
          );
        })}
      </tbody>
    </table>
  );
}
