"""make reserved submission codes unique

Revision ID: b731f4c8d2a1
Revises: e38e216beb58
"""

from alembic import op


revision = "b731f4c8d2a1"
down_revision = "e38e216beb58"
branch_labels = None
depends_on = None


def upgrade():
    op.create_index(
        "uq_submission_proposed_code",
        "submission",
        ["proposed_code"],
        unique=True,
    )


def downgrade():
    op.drop_index(
        "uq_submission_proposed_code",
        table_name="submission",
    )
