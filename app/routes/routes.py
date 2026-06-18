import os
import tempfile
import datetime
import base64
import json

from werkzeug.utils import secure_filename
from app.chemistry.processor import process_structure

# importing important route-related flask functions, form for searching database, database models, blueprint for routes
from flask import render_template, redirect, url_for, abort, flash, make_response, send_file
from flask_login import current_user, login_user, logout_user, login_required
from app.routes.forms import SearchForm, ImportPeptoidForm, ContributorLoginForm, SubmitSubmissionForm
from app.models import Peptoid, Author, Residue, Contributor, Submission
from app.routes import bp
from app import app, db
from flask import request

# function for creating all gallery views


def renderGallery(peptoids, title, page, next_url, prev_url, view, var):
    if view not in ['3d', '2d']:
        flash(f'INVALID VIEW ARGUMENT: {view}. Must be \'3d\' or \'2d\'','danger')
        abort(400)
    
    sequence_max = 128  # sequence cap

    # lists passed to generate gallery
    peptoid_codes = []
    peptoid_urls = []
    peptoid_titles = []
    peptoid_sequences = []
    data = []
    publications = []

    # populating lists with properties from peptoids passed by route function
    for p in peptoids.items:
        peptoid_codes.append(p.code)
        peptoid_titles.append(p.title)
        peptoid_urls.append(url_for('routes.peptoid', code=p.code))

        # creating peptoid sequence string according adhering to the max sequence characters (either full sequence or first 3 residues and ...)
        seq_val = p.sequence or ""
        if len(seq_val) < sequence_max:
            peptoid_sequences.append(seq_val)
        else:
            i = [pos for pos, char in enumerate(seq_val) if char == ',']
            if len(i) > 2:
                peptoid_sequences.append(seq_val[:i[2]] + " ...")
            else:
                peptoid_sequences.append(seq_val[:sequence_max] + " ...")

        # if structure doi exists add doi links for structure doi and pub doi, else use empty string for struct doi link
        if p.struct_doi:
            data.append("https://www.doi.org/{}".format(p.struct_doi))
            publications.append("https://www.doi.org/{}".format(p.pub_doi))
        else:
            data.append("")
            publications.append("https://www.doi.org/{}".format(p.pub_doi))

    # return gallery.html template with properties
    return render_template('gallery.html',
                           title=title,
                           pages=peptoids.pages,
                           total=peptoids.total,
                           next_url=next_url,
                           prev_url=prev_url,
                           page=page,
                           view=view,
                           var=var,
                           peptoid_codes=peptoid_codes,
                           peptoid_urls=peptoid_urls,
                           peptoid_titles=peptoid_titles,
                           peptoid_sequences=peptoid_sequences,
                           data=data,
                           publications=publications
                           )


@bp.route('/')
@bp.route('/home')
def home():
    return render_template('home.html', title="Peptoid Data Bank")

# route for baseline gallery, passes var=False because no var argument


@bp.route('/gallery')
def gallery():
    title = 'Gallery'
    # page number argument for pagination
    page = request.args.get('page', 1, type=int)
    # view argument for 3d vs 2d
    view = request.args.get('view', '2d', type=str)

    # paginated peptoids
    peptoids = Peptoid.query.order_by(Peptoid.release.desc()).paginate(
    page=page, per_page=app.config['PEPTOIDS_PER_PAGE'], error_out=False)
    next_url = url_for(
        'routes.gallery', page=peptoids.next_num, view=view) if peptoids.has_next else None
    prev_url = url_for(
        'routes.gallery', page=peptoids.prev_num, view=view) if peptoids.has_prev else None

    return renderGallery(peptoids, title, page, next_url, prev_url, view, var=False)

# route for all residues renders residues.html


@bp.route('/residues')
def residues():
    title = 'Residues'
    residues = [r for r in Residue.query.all()]
    return render_template('residues.html', title=title, residues=residues)

# search route renders the form using the search.html template and the SearchForm() from forms.py
# if the user has submitted the form they are redirected to the route for the serial choice with
# var = their search box input


@bp.route('/search', methods=['GET', 'POST'])
def search():
    form = SearchForm()

    if request.method == 'POST':
        filters = {}

        for key in ['residue', 'author', 'doi', 'experiment', 'topology']:
            value = getattr(form, key).data
            if value:
                value = value.strip() if isinstance(value, str) else value
                if value:
                    filters[key] = value

        if not filters:
            flash('Please enter at least one search value.', 'warning')
            return render_template(
                'search.html',
                title='Search',
                form=form,
                info='Filter database of peptoids by one or more of the properties below.',
                description='Peptoid Data Bank - Explore by property'
            )

        for key, value in list(filters.items()):
            if key == 'doi':
                filters[key] = value.replace('/', '$')

        return redirect(url_for('routes.multisearch', **filters))

    return render_template(
        'search.html',
        title='Search',
        form=form,
        info='Filter database of peptoids by one or more of the properties below.',
        description='Peptoid Data Bank - Explore by property'
    )


@bp.route('/multisearch')
def multisearch():
    page = request.args.get('page', 1, type=int)
    view = request.args.get('view', '2d', type=str)

    filters = {}
    for key in ['residue', 'author', 'doi', 'experiment', 'topology']:
        value = request.args.get(key, type=str)
        if value:
            filters[key] = value

    if not filters:
        flash('No search filters were provided.', 'warning')
        return redirect(url_for('routes.search'))

    matching_codes = []

    for p in Peptoid.query.all():
        match = True

        for key, raw_value in filters.items():
            value = raw_value.replace('$', '/').strip().lower()

            if key == 'residue':
                residue_values = []
                for r in p.peptoid_residue:
                    residue_values.append((r.short_name or '').lower())
                    residue_values.append((r.long_name or '').lower())

                if not any(value in item for item in residue_values):
                    match = False
                    break

            elif key == 'author':
                author_values = []
                for a in p.peptoid_author:
                    first = (a.first_name or '').lower()
                    last = (a.last_name or '').lower()
                    author_values.extend([
                        first,
                        last,
                        f"{first} {last}".strip(),
                        f"{last}, {first}".strip()
                    ])

                if not any(value in item for item in author_values):
                    match = False
                    break

            elif key == 'doi':
                doi_values = [
                    (p.struct_doi or '').lower(),
                    (p.pub_doi or '').lower()
                ]

                if not any(value in item for item in doi_values):
                    match = False
                    break

            elif key == 'experiment':
                if value != (p.experiment or '').lower():
                    match = False
                    break

            elif key == 'topology':
                if value != (p.topology or '').lower():
                    match = False
                    break

        if match:
            matching_codes.append(p.code)

    peptoids = Peptoid.query.filter(Peptoid.code.in_(matching_codes)).order_by(
        Peptoid.release.desc()
    ).paginate(
        page=page,
        per_page=app.config['PEPTOIDS_PER_PAGE'],
        error_out=False
    )

    next_url = url_for(
        'routes.multisearch',
        page=peptoids.next_num,
        view=view,
        **filters
    ) if peptoids.has_next else None

    prev_url = url_for(
        'routes.multisearch',
        page=peptoids.prev_num,
        view=view,
        **filters
    ) if peptoids.has_prev else None

    labels = {
        'residue': 'Residue',
        'author': 'Author',
        'doi': 'DOI',
        'experiment': 'Experiment',
        'topology': 'Topology'
    }

    display_filters = []
    for key, value in filters.items():
        shown_value = value.replace('$', '/')
        display_filters.append(f"{labels.get(key, key.title())}: {shown_value}")

    title = 'Search Results: ' + ', '.join(display_filters)

    return renderGallery(peptoids, title, page, next_url, prev_url, view, filters)


# individual peptoid page for peptoid specified by code

# individual peptoid page for peptoid specified by code


@bp.route('/peptoid/<code>')
def peptoid(code):
    # data passed to front end
    peptoid = Peptoid.query.filter_by(code=code).first_or_404()
    title = peptoid.title
    code = peptoid.code
    release = str(peptoid.release.month) + "/" + \
        str(peptoid.release.day) + "/" + str(peptoid.release.year)
    experiment = peptoid.experiment

    # lists of objects also passed to front end
    authors = []
    residues = []

    for author in peptoid.peptoid_author:
        authors.append(author)

    for residue in peptoid.peptoid_residue:
        residues.append(residue)

    sequence = peptoid.sequence
    author_list = ", ".join(
        [a.first_name + " " + a.last_name for a in authors])

    data = peptoid.struct_doi
    publication = peptoid.pub_doi

    return render_template('peptoid.html',
                           peptoid=peptoid,
                           title=title,
                           code=code,
                           release=release,
                           experiment=experiment,
                           data=data,
                           publication=publication,
                           authors=authors,
                           residues=residues,
                           sequence=sequence,
                           author_list=author_list
                           )

# popout for residues of peptoid


@bp.route('/residue/<var>/popout')
def residue_popout(var):
    residue = Residue.query.filter((Residue.long_name == var) | (
        Residue.short_name == var)).first_or_404()
    return render_template('residue_popout.html', residue=residue, title='Popout')

@bp.route('/peptoid/<code>/citation')
def citation(code):
    peptoid = Peptoid.query.filter_by(code=code).first_or_404()
    response = make_response(peptoid.citation, 200)
    response.mimetype = "text/plain"
    return response


# residue route for residue, name = var
@bp.route('/residue/<var>')
def residue(var):
    page = request.args.get('page', 1, type=int)
    view = request.args.get('view', '2d', type=str)
    title = 'Filtered by Residue: ' + var
    # get residue(s) and related peptoids codes
    residues = Residue.query.filter((Residue.long_name == var) | (
            Residue.short_name == var)).all()
    if len(residues) == 0:
        flash(f'NO RESIDUES for: <{var}>','danger')
        abort(404)
    peptoids = []
    for r in residues:
        peptoids.extend(r.peptoids)
    codes = [p.code for p in peptoids]
    peptoids = Peptoid.query.filter(Peptoid.code.in_(codes)).order_by(Peptoid.release.desc()).paginate(
        page, app.config['PEPTOIDS_PER_PAGE'], True)  # querying for peptoids based on list of codes
    next_url = url_for('routes.residue', page=peptoids.next_num,
                       var=var, view=view) if peptoids.has_next else None
    prev_url = url_for('routes.residue', page=peptoids.prev_num,
                       var=var, view=view) if peptoids.has_prev else None
    return renderGallery(peptoids, title, page, next_url, prev_url, view, var)

# author route for author. If name entered has a space search by both words for first name and last name.
# if name entered is just one word check if it is an author's first name or last name


@bp.route('/author/<var>')
def author(var):
    page = request.args.get('page', 1, type=int)
    view = request.args.get('view', '2d', type=str)

    # retrieving peptoids from author according to input var, processed for comma separated name
    if "," in var:
        name_split = var.split(', ')
        last_name = name_split[0]
        first_name = name_split[1]
        authors = Author.query.filter_by(
            first_name=first_name, last_name=last_name).all()
    else:
        authors = Author.query.filter((Author.first_name == var) | (
            Author.last_name == var)).all()
    if len(authors) == 0:
        flash(f'NO AUTHORS for: <{var}>','danger')
        abort(404)
    peptoids = []
    for a in authors:
        peptoids.extend(a.peptoids)
    codes = [p.code for p in peptoids]
    title = 'Filtered by Author: ' + var
    peptoids = Peptoid.query.filter(Peptoid.code.in_(codes)).order_by(Peptoid.release.desc()).paginate(
        page, app.config['PEPTOIDS_PER_PAGE'], True)
    next_url = url_for('routes.author', page=peptoids.next_num,
                       var=var, view=view) if peptoids.has_next else None
    prev_url = url_for('routes.author', page=peptoids.prev_num,
                       var=var, view=view) if peptoids.has_prev else None
    return renderGallery(peptoids, title, page, next_url, prev_url, view, var)

# experiment route for Peptoid.experimet = var, returns gallery.html (gallery view), if no peptoid found returns a 404


@bp.route('/experiment/<var>')
def experiment(var):
    page = request.args.get('page', 1, type=int)
    view = request.args.get('view', '2d', type=str)
    title = 'Filtered by Experiment: ' + var
    peptoids = Peptoid.query.order_by(Peptoid.release.desc()).filter_by(
        experiment=var).paginate(page, app.config['PEPTOIDS_PER_PAGE'], True)
    if len(peptoids.items) == 0:
        abort(404)
    next_url = url_for('routes.experiment', page=peptoids.next_num,
                       var=var, view=view) if peptoids.has_next else None
    prev_url = url_for('routes.experiment', page=peptoids.prev_num,
                       var=var, view=view) if peptoids.has_prev else None
    return renderGallery(peptoids, title, page, next_url, prev_url, view, var)

# doi route for Peptoid.doi = var, returns gallery.html (gallery view), if no peptoid found returns a 404


@bp.route('/doi/<var>')
def doi(var):
    page = request.args.get('page', 1, type=int)
    var = var.replace('$', '/') #making doi into doi that fits database
    view = request.args.get('view', '2d', type=str)
    title = 'Filtered by DOI: ' + var
    peptoids = Peptoid.query.order_by(Peptoid.release.desc()).filter(
        (Peptoid.struct_doi == var) | (Peptoid.pub_doi == var)).paginate(page, app.config['PEPTOIDS_PER_PAGE'], True)
    if len(peptoids.items) == 0:
        abort(404)
    var = var.replace('/', '$') #making doi fit url
    next_url = url_for('routes.doi', page=peptoids.next_num,
                       var=var, view=view) if peptoids.has_next else None
    prev_url = url_for('routes.doi', page=peptoids.prev_num,
                       var=var, view=view) if peptoids.has_prev else None
    return renderGallery(peptoids, title, page, next_url, prev_url, view, var)

# filtering according to topology


@bp.route('/top/<var>')
def topology(var):
    page = request.args.get('page', 1, type=int)
    view = request.args.get('view', '2d', type=str)
    title = 'Filtered by Topology: ' + var
    peptoids = Peptoid.query.order_by(Peptoid.release.desc()).filter_by(
        topology=var).paginate(page, app.config['PEPTOIDS_PER_PAGE'], True)
    if len(peptoids.items) == 0:
        abort(404)
    next_url = url_for('routes.topology', page=peptoids.next_num,
                       var=var, view=view) if peptoids.has_next else None
    prev_url = url_for('routes.topology', page=peptoids.prev_num,
                       var=var, view=view) if peptoids.has_prev else None
    return renderGallery(peptoids, title, page, next_url, prev_url, view, var)


# api route returns api.html template
@bp.route('/api', methods=['GET', 'POST'])
def api():
    return render_template('api.html', title="PeptoidDB API")

@bp.route('/contribute')
def contribute():
    submissions = []
    if current_user.is_authenticated:
        submissions = Submission.query.filter_by(
            contributor_id=current_user.id
        ).order_by(Submission.updated_at.desc()).all()

    return render_template(
        'contribute.html',
        title='Contribute',
        submissions=submissions,
    )


@bp.route('/submission/<int:submission_id>')
@login_required
def view_submission(submission_id):
    submission = Submission.query.filter_by(
        id=submission_id,
        contributor_id=current_user.id,
    ).first_or_404()

    return render_template(
        'submission.html',
        title='Submission {}'.format(submission.id),
        submission=submission,
        residues=json.loads(submission.residue_data_json),
        submit_form=SubmitSubmissionForm(),
    )


@bp.route('/submission/<int:submission_id>/image/<kind>')
@login_required
def submission_image(submission_id, kind):
    submission = Submission.query.filter_by(
        id=submission_id,
        contributor_id=current_user.id,
    ).first_or_404()

    paths = {
        'structure': submission.structure_image_path,
        'residues': submission.residue_image_path,
    }
    path = paths.get(kind)

    if not path or not os.path.isfile(path):
        abort(404)

    return send_file(path, mimetype='image/png')


@bp.route('/contributor-login', methods=['GET', 'POST'])
def contributor_login():
    if current_user.is_authenticated:
        return redirect(url_for('routes.contribute'))

    form = ContributorLoginForm()
    if form.validate_on_submit():
        contributor = Contributor.query.filter_by(username=form.username.data.strip()).first()
        if contributor and contributor.is_active and contributor.check_password(form.password.data):
            login_user(contributor, remember=form.remember.data)
            contributor.last_login = datetime.datetime.utcnow()
            db.session.commit()
            next_page = request.args.get('next')
            if not next_page or not next_page.startswith('/') or next_page.startswith('//'):
                next_page = url_for('routes.contribute')
            return redirect(next_page)
        flash('Invalid username or password.', 'danger')

    return render_template('contributor_login.html', title='Contributor Login', form=form)


@bp.route('/contributor-logout')
@login_required
def contributor_logout():
    logout_user()
    flash('You have been logged out.', 'info')
    return redirect(url_for('routes.contribute'))


@bp.route('/import-peptoid', methods=['GET', 'POST'])
@login_required
def import_peptoid():
    form = ImportPeptoidForm()
    submit_form = SubmitSubmissionForm()
    preview = None
    submission = None

    if form.validate_on_submit():
        try:
            cif_bytes = None
            original_cif_filename = None

            if form.cif_file.data and form.cif_file.data.filename:
                original_cif_filename = secure_filename(
                    form.cif_file.data.filename
                )
                cif_bytes = form.cif_file.data.read()

                with tempfile.TemporaryDirectory() as temp_dir:
                    cif_path = os.path.join(
                        temp_dir, original_cif_filename
                    )
                    with open(cif_path, 'wb') as cif_file:
                        cif_file.write(cif_bytes)
                    preview = process_structure(cif_path=cif_path)
            else:
                preview = process_structure(smiles=form.smiles.data)

            release_datetime = datetime.datetime.combine(
                form.release.data, datetime.time.min
            )
            submission = Submission(
                contributor_id=current_user.id,
                status='draft',
                proposed_code=form.code.data.strip(),
                title=form.title.data.strip(),
                release=release_datetime,
                experiment=form.experiment.data,
                pub_doi=(form.pub_doi.data or '').strip() or None,
                struct_doi=(form.struct_doi.data or '').strip() or None,
                citation=(form.citation.data or '').strip() or None,
                authors=form.authors.data.strip(),
                input_type=preview['input_type'],
                original_smiles=(
                    form.smiles.data.strip()
                    if preview['input_type'] == 'SMILES' else None
                ),
                cleaned_smiles=preview['cleaned_smiles'],
                topology=preview['topology'],
                residue_data_json=json.dumps(preview['residues']),
                warnings_json=json.dumps([]),
                original_cif_filename=original_cif_filename,
            )
            db.session.add(submission)
            db.session.flush()

            staging_dir = os.path.join(
                app.instance_path,
                'submission_staging',
                str(submission.id),
            )
            os.makedirs(staging_dir, exist_ok=True)

            structure_path = os.path.join(staging_dir, 'structure.png')
            residue_path = os.path.join(staging_dir, 'residues.png')

            with open(structure_path, 'wb') as image_file:
                image_file.write(
                    base64.b64decode(preview['structure_image'])
                )
            with open(residue_path, 'wb') as image_file:
                image_file.write(
                    base64.b64decode(preview['residue_image'])
                )

            submission.structure_image_path = structure_path
            submission.residue_image_path = residue_path

            if cif_bytes is not None:
                staged_cif_path = os.path.join(
                    staging_dir, original_cif_filename
                )
                with open(staged_cif_path, 'wb') as cif_file:
                    cif_file.write(cif_bytes)
                submission.staged_cif_path = staged_cif_path

            db.session.commit()
            flash(
                'Draft saved. Review the preview, then submit it for review.',
                'success',
            )

        except Exception as error:
            db.session.rollback()
            flash(f'Could not process structure: {error}', 'danger')

    return render_template(
        'import_peptoid.html',
        title='Import Peptoid',
        form=form,
        submit_form=submit_form,
        preview=preview,
        submission=submission,
    )


@bp.route('/submission/<int:submission_id>/submit', methods=['POST'])
@login_required
def submit_submission(submission_id):
    form = SubmitSubmissionForm()
    submission = Submission.query.filter_by(
        id=submission_id,
        contributor_id=current_user.id,
    ).first_or_404()

    if not form.validate_on_submit():
        abort(400)

    if submission.status != 'draft':
        flash('This submission is no longer a draft.', 'warning')
        return redirect(url_for('routes.contribute'))

    submission.status = 'pending'
    db.session.commit()
    flash('Your entry was submitted for review.', 'success')
    return redirect(url_for('routes.contribute'))
